#include "compute_lk.h"
#include <mpi.h>
#include <cstring>
#include "fix_store.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "error.h"
#include "memory.h"
#include "domain.h"
#include "group.h"
#include "math_extra.h"
#include "math_const.h"
#include "compute_fragment_atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace LAMMPS_NS;

ComputeLK::ComputeLK(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), list(NULL)
{
  if (narg > 10) error->all(FLERR, "Illegal compute lk command");   //parse arguments
  for (int iarg=3; iarg<narg; ++iarg)
  {
    if (strcmp(arg[iarg], "nonperiodic") == 0)
    {
      pflag=false;
    }
    else if (strcmp(arg[iarg], "periodic") == 0)
    {
      pflag=true;
    }
    else if (strcmp(arg[iarg], "molecular") == 0)
    {
      frflag=false;
    }
    else if (strcmp(arg[iarg], "fragment") == 0)    //find compute/fragment command
    {
      frflag=true;
      int p=strlen(arg[iarg+1])+1;
      idfragment=new char[p];
      strcpy(idfragment, arg[iarg+1]);
      int icompute=modify->find_compute(idfragment);
      if (icompute<0)
      {
        error->all(FLERR, "fragment/atom compute or group does not exist for compute lk");
      }
      cfragment=(ComputeFragmentAtom *) modify->compute[icompute];
      if (strcmp(cfragment->style, "fragment/atom") != 0)
      {
        error->all(FLERR, "compute lk does not use fragment/atom compute or group");
      }
      ++iarg;
    }
    else if (strcmp(arg[iarg], "writhe") == 0)
    {
      wrflag=true;
    }
    else if (strcmp(arg[iarg], "linking") == 0)
    {
      wrflag=false;
    }
    else if (strcmp(arg[iarg], "group2") == 0)    //find group2
    {
      int m=strlen(arg[iarg+1])+1;
      group2=new char[m];
      strcpy(group2,arg[iarg+1]);
      jgroup=group->find(group2);
      if (jgroup > -1)
      {
        jgroupbit=group->bitmask[jgroup];
        grflag=true;
      }
      else
      {
        error->all(FLERR, "Illegal group2");
      }
      ++iarg;
    }
    else
    {
      error->all(FLERR, "Illegal argument");
    }
    
  }

  /*int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;    //save positions of atoms
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";
  newarg[5] = (char *) "3";
  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;

    double **x = atom->x;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
        xoriginal[i][0]=x[i][0];
        xoriginal[i][1]=x[i][1];
        xoriginal[i][2]=x[i][2];
    }
  }*/

  if (wrflag)
  {
    vector_flag=1;
    extvector=0;
    array_flag=0;
  }
  else
  {
    array_flag=1;
    extarray=0;
    vector_flag=0;
  }
  

  tagint *molecule=atom->molecule;
  tagint *tag = atom->tag;
  int nlocal=atom->nlocal;
  int *mask=atom->mask;
  mols=NULL;
  nmols=0;
  int atom1id, atom2id;
  if (!frflag) {       //count number of chains if molecular
    for (int i=0; i<nlocal; ++i)
    {
      atom1id = tag[i]-1;
      found=false;
      if (!(mask[atom1id] & groupbit) && !(mask[atom1id] & jgroupbit) && grflag) continue;
      if (!frflag) {
        for (int j=0; j<nmols; ++j)
        {
          atom2id = tag[j]-1;
          if (mols[atom2id] == molecule[atom1id]) found=true;
        }
        if (!found)
        {
          memory->grow(mols, nmols+1, "mols");
          mols[nmols]=molecule[atom1id];
          nmols++;
        }
      }
    }
  }
  else {      //count number of chains if fragment
    cfragment->compute_peratom();
    double *fragment_vec=cfragment->vector_atom;
    memory->create(frmolecule, nlocal, "frmolecule");
    int ct=0;
    for (int i=0; i<nlocal; ++i)
    {
      atom1id = tag[i]-1;
      if (ct == 0 || fragment_vec[atom1id] != fragment_vec[atom1id-1])
      {
        ++ct;
        memory->grow(mols, ct, "mols");
        mols[ct-1]=ct;
      }
      frmolecule[atom1id]=(tagint) ct;
    }
    nmols=ct;
  }
  if (wrflag)
  {
    size_vector=nmols;
    vector=new double[nmols];
  }
  else
  {
    size_array_cols=nmols;
    size_array_rows=nmols;
    memory->create(array,size_array_rows,size_array_cols,"lk:array");
  }
}

ComputeLK::~ComputeLK()
{
  //if (modify->nfix) modify->delete_fix(id_fix);
  if (array_flag==1) memory->destroy(array);
  if (vector_flag==1) 
  {
    delete [] vector;
  }
  if (frflag)
  {
    delete [] idfragment;
  }
  delete [] group2;
}

void ComputeLK::allocate()  //NOT USED
{/*
  tagint *molecule=atom->molecule;
  int nlocal=atom->nlocal;
  int *mask=atom->mask;
  mols=NULL;
  nmols=0;
  if (!frflag) {
    for (int i=0; i<nlocal; ++i)
    {
      found=false;
      if (!(mask[i] & groupbit) && !(mask[i] & jgroupbit) && grflag) continue;
      if (!frflag) {
        for (int j=0; j<nmols; ++j)
        {
          if (mols[j] == molecule[i]) found=true;
        }
        if (!found)
        {
          memory->grow(mols, nmols+1, "mols");
          mols[nmols]=molecule[i];
          nmols++;
        }
      }
    }
  }
  else {
    cfragment->compute_peratom();
    double *fragment_vec=cfragment->vector_atom;
    memory->grow(molecule, nlocal, "molecule");
    int ct=0;
    for (int i=0; i<nlocal; ++i)
    {
      if (ct == 0 || fragment_vec[i] != fragment_vec[i-1])
      {
        ++ct;
        memory->grow(mols, ct, "mols");
        mols[ct-1]=ct;
      }
      molecule[i]=(tagint) ct;
    }
    nmols=ct;
  }

  if (!wrflag)
  {
    size_array_cols=size_array_rows=nmols;
    memory->create(array,size_array_rows,size_array_cols,"lk:array");
  }*/
}

void ComputeLK::init()
{
  /*int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute msd fix ID");
  fix = (FixStore *) modify->fix[ifix];*/
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

void ComputeLK::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

void ComputeLK::compute_vector()
{
  invoked_vector=update->ntimestep;
  double **x=atom->x;
  tagint *molecule=atom->molecule;
  tagint *tag = atom->tag;
  int nlocal=atom->nlocal;
  int nmax=atom->nmax;
  int *mask=atom->mask;
  imageint imagecpy[nlocal];
  imageint *image=atom->image;
  int xbox, ybox, zbox;
  int trans[3];
  int tmpx,tmpy,tmpz;
  imageint tmp;
  double lkp=0;
  int ch1=0;
  int ch2=0;
  double results;
  if (frflag)
  {
    memory->grow(molecule, nlocal, "molecule");
    for (int i=0; i<nlocal; ++i)
    {
      molecule[i]=frmolecule[i];
    }
  }

  if (pflag)
  {
    for (int i = 0; i < nlocal; i++)
    {
      imagecpy[tag[i]-1] = image[i];
    }
  }

  double dx, dy, dz;
  int closest;
  int atomid;

  for (int k=0; k<nmols; k++)
  {
    ch1=mols[k];
    imageint *transimg=NULL;
    imageint *img1=NULL;
    imageint *img2=NULL;
    int ntrans=0;
    int nimg1=0;
    int nimg2=0;
    ch2=ch1;
    lkp=0;
    int ct=0;
    int atom1id = 0;
    int atom2id = 0;
    for (int i=0; i<nlocal; ++i)
    {
      atom1id = tag[i]-1;
      if (molecule[i] != ch1) continue;
    
      if (pflag)    //create array images the chain touches
      {
        found=false;

        for (int j=0; j<nimg1; ++j)
        {
          if (img1[j]==imagecpy[i]) found=true;
        }
        if (!found)
        {
          memory->grow(img1, nimg1+1, "img1");
          img1[nimg1]=imagecpy[i];
          nimg1++;
        }
      }
    }
    if (pflag)
    {
      int ixbox,iybox,izbox,jxbox,jybox,jzbox;
      for (int i=0; i<nimg1; ++i)               //extract image details from packed bits
      {
        ixbox = (img1[0] & IMGMASK) - IMGMAX;
        iybox = (img1[0] >> IMGBITS & IMGMASK) - IMGMAX;
        izbox = (img1[0] >> IMG2BITS) - IMGMAX;
          
        jxbox = (img1[i] & IMGMASK) - IMGMAX;
        jybox = (img1[i] >> IMGBITS & IMGMASK) - IMGMAX;
        jzbox = (img1[i] >> IMG2BITS) - IMGMAX;
    
        tmp = ((imageint) ((ixbox-jxbox) + IMGMAX) & IMGMASK) |           //subtract image2 from image 1 and convert back to packed bits
          (((imageint) ((iybox-jybox) + IMGMAX) & IMGMASK) << IMGBITS) |
          (((imageint) ((izbox-jzbox) + IMGMAX) & IMGMASK) << IMG2BITS);
    
        found=false;
        for (int j=0; j<ntrans; ++j)
        {
          if (transimg[j]==tmp) found=true;
        }
    
        if (!found)   //create array of unique translated images and compute the lk of that translation
        {
          memory->grow(transimg, ntrans+1, "transimg");
          transimg[ntrans]=tmp;
          ntrans++;
          lkp+=lk(ch1,ch2,(ixbox-jxbox)*domain->xprd, (iybox-jybox)*domain->yprd, (izbox-jzbox)*domain->zprd);
        }  
      }
      vector[k]=lkp;
    }
      
    else    //calculate nonperiodic writhe if desired
    {
      vector[k]=lk(ch1,ch2);
    }
  }
}

void ComputeLK::compute_array()   //same as above except for lk instead of writhe
{
  invoked_array=update->ntimestep;
  double **x=atom->x;
  tagint *molecule=atom->molecule;
  tagint *tag = atom->tag;
  int nlocal=atom->nlocal;
  int nmax=atom->nmax;
  int *mask=atom->mask;
  imageint imagecpy[nlocal];
  imageint *image=atom->image;
  int xbox, ybox, zbox;
  int trans[3];
  int tmpx,tmpy,tmpz;
  imageint tmp;
  double lkp=0;
  int ch1=0;
  int ch2=0;
  double temp[2][2];

  neighbor->build_one(list);
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  std::cout << "inum: " << inum << "\n";
  for (int i = 0; i < inum; i++) {
    std::cout << "ilist[" << i << "] = " << ilist[i] << "\n";
  }
  if (frflag)
  {
    memory->grow(molecule, nlocal, "molecule");
    for (int i=0; i<nlocal; ++i)
    {
      molecule[i]=frmolecule[i];
    }
  }

  if (pflag)
  {
    for (int i = 0; i < nlocal; i++)
    {
      imagecpy[tag[i]-1] = image[i];
    }
  }

  double dx, dy, dz;
  int closest;
  int atomid;
  maximum = 0;

  for (int k=0; k<nmols; k++)
  {
    ch1=mols[k];
    imageint *transimg=NULL;
    imageint *img1=NULL;
    imageint *img2=NULL;
    int ntrans=0;
    int nimg1=0;
    int nimg2=0;
    ch2=ch1;
    lkp=0;
    int ct=0;
    int atom1id = 0;
    int atom2id = 0;
    for (int i=0; i<nlocal; ++i)
    {
      atom1id = tag[i]-1;
      if (molecule[i] != ch1) continue;
    
      if (pflag)    //create array images the chain touches
      {
        found=false;

        for (int j=0; j<nimg1; ++j)
        {
          if (img1[j]==imagecpy[i]) found=true;
        }
        if (!found)
        {
          memory->grow(img1, nimg1+1, "img1");
          img1[nimg1]=imagecpy[i];
          int ixbox = (imagecpy[i] & IMGMASK) - IMGMAX;
          int iybox = (imagecpy[i] >> IMGBITS & IMGMASK) - IMGMAX;
          int izbox = (imagecpy[i] >> IMG2BITS) - IMGMAX;
          nimg1++;
        }
      }
    }
    if (pflag)
    {
      int ixbox,iybox,izbox,jxbox,jybox,jzbox;
      for (int i=0; i<nimg1; ++i)               //extract image details from packed bits
      {
        ixbox = (img1[0] & IMGMASK) - IMGMAX;
        iybox = (img1[0] >> IMGBITS & IMGMASK) - IMGMAX;
        izbox = (img1[0] >> IMG2BITS) - IMGMAX;
          
        jxbox = (img1[i] & IMGMASK) - IMGMAX;
        jybox = (img1[i] >> IMGBITS & IMGMASK) - IMGMAX;
        jzbox = (img1[i] >> IMG2BITS) - IMGMAX;
    
        tmp = ((imageint) ((ixbox-jxbox) + IMGMAX) & IMGMASK) |           //subtract image2 from image 1 and convert back to packed bits
          (((imageint) ((iybox-jybox) + IMGMAX) & IMGMASK) << IMGBITS) |
          (((imageint) ((izbox-jzbox) + IMGMAX) & IMGMASK) << IMG2BITS);
    
        found=false;
        for (int j=0; j<ntrans; ++j)
        {
          if (transimg[j]==tmp) found=true;
        }
    
        if (!found)   //create array of unique translated images and compute the lk of that translation
        {
          memory->grow(transimg, ntrans+1, "transimg");
          transimg[ntrans]=tmp;
          ntrans++;
          lkp+=lk(ch1,ch2,(ixbox-jxbox)*domain->xprd, (iybox-jybox)*domain->yprd, (izbox-jzbox)*domain->zprd);
          //std::cout << k << " " << k << " " << lkp << std::endl;
        }  
      }
      array[k][k]=lkp;
    }
      
    else    //calculate nonperiodic writhe if desired
    {
      array[k][k]=lk(ch1,ch2);
    }

    for (int m=k+1; m<nmols; m++)
    {
      ch2=mols[m];
      transimg=NULL;
      img1=NULL;
      img2=NULL;
      ntrans=0;
      nimg1=0;
      nimg2=0;
      lkp=0;
      int ct=0;

      for (int i=0; i<nlocal; ++i)
      {
        if (molecule[i] != ch1 && molecule[i] != ch2) continue;
      
        if (pflag)    //generate arrays of images that each chain touch
        {
          atom1id = tag[i]-1;
          found=false;
          if (molecule[i]==ch1)
          {
            for (int j=0; j<nimg1; ++j)
            {
              if (img1[j]==imagecpy[i]) found=true;
            }
            if (!found)
            {
              memory->grow(img1, nimg1+1, "img1");
              img1[nimg1]=imagecpy[i];
              nimg1++;
            }
          }
          else
          {
            for (int j=0; j<nimg2; ++j)
            {
              if (img2[j]==imagecpy[i]) found=true;
            }
            if (!found)
            {
              memory->grow(img2, nimg2+1, "img2");
              img2[nimg2]=imagecpy[i];
              nimg2++;
            }
          }
        }
      }
          
      if (pflag)
      {
        int ixbox,iybox,izbox,jxbox,jybox,jzbox;
        for (int i=0; i<nimg1; ++i)
        {
          ixbox = (img1[i] & IMGMASK) - IMGMAX;
          iybox = (img1[i] >> IMGBITS & IMGMASK) - IMGMAX;
          izbox = (img1[i] >> IMG2BITS) - IMGMAX;
          for (int j=0; j<nimg2; j++)
          {
            jxbox = (img2[j] & IMGMASK) - IMGMAX;
            jybox = (img2[j] >> IMGBITS & IMGMASK) - IMGMAX;
            jzbox = (img2[j] >> IMG2BITS) - IMGMAX;
            /*if (ch1 == 1 && ch2 == 2)
            {
              std::cout << ixbox << " " << iybox << " " << izbox << std::endl;
              std::cout << jxbox << " " << jybox << " " << jzbox << std::endl;
              std::cout << " " << std::endl;
              std::cout << ixbox-jxbox << " " << iybox-jybox << " " << izbox-jzbox << std::endl;
            }*/
      
            tmp = ((imageint) ((ixbox-jxbox) + IMGMAX) & IMGMASK) |
              (((imageint) ((iybox-jybox) + IMGMAX) & IMGMASK) << IMGBITS) |
              (((imageint) ((izbox-jzbox) + IMGMAX) & IMGMASK) << IMG2BITS);
      
            found=false;
            for (int l=0; l<ntrans; ++l)
            {
              if (transimg[l]==tmp) found=true;
            }
      
            if (!found)
            {
              memory->grow(transimg, ntrans+1, "transimg");   //create array of unique translations of the chains
              transimg[ntrans]=tmp;
              ntrans++;
              double temp_lk = lk(ch1,ch2,(ixbox-jxbox)*domain->xprd, (iybox-jybox)*domain->yprd, (izbox-jzbox)*domain->zprd);
              lkp+=temp_lk;
              if (abs(temp_lk) > maximum) {
                  maximum = abs(temp_lk);
                  maxch1 = ch1;
                  maxch2 = ch2;
                  maxxtrans = (ixbox-jxbox)*domain->xprd;
                  maxytrans = (iybox-jybox)*domain->yprd;
                  maxztrans = (izbox-jzbox)*domain->zprd;
              }
              ct++;
            }
          }
        }
        array[k][m]=array[m][k]=lkp;
      }
        
      else 
      {
        array[k][m]=lk(ch1,ch2);    //calculate nonperiodic lk if desired
        array[m][k]=array[k][m];
      } 
    }
  }

  std::ofstream avgfile;
  avgfile.open("avgfile.txt");
  double totallk=0;
  double totalwr=0;
  double lkiter=0;
  double writer=0;
  for (int i=0; i<nmols; ++i)
  {
    for (int j=0; j<nmols; ++j)
    {
      if (i != j) {
        totallk+=abs(array[i][j]);
        //std::cout << abs(array[i][j]) << std::endl;
        ++lkiter;
      }

      else {
        totalwr+=abs(array[i][j]);
        ++writer;
      }
    }
  }
  double avlk = totallk/lkiter;
  double avwr = totalwr/writer;
  avgfile << "Average linking number (absolute value): " << avlk << "\n";
  avgfile << "Average writhe (absolute value): " << avwr << "\n";
  avgfile << "Maximum absolute linking number: " << maximum << "\n";
  avgfile << "Chain1: " << maxch1 << ", chain2: " << maxch2 << "\n";
  avgfile << "xtrans: " << maxxtrans << ", ytrans: " << maxytrans << ", ztrans: " << maxztrans << "\n";
  avgfile.close();
}

double ComputeLK::lk(int ch1, int ch2, double jxtrans, double jytrans, double jztrans)    //method computes lk between two given chains
{
  tagint *molecule=atom->molecule;    
  int nlocal=atom->nlocal;
  double **x=atom->x;
  int **bond_type=atom->bond_type;
  tagint **bond_atom=atom->bond_atom;
  imageint *image=atom->image;
  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  //int i1, i2, j1, j2, atom1id, atom2id;
  int f1=-1;
  int f2=-1;
  double lk=0;
  double xcpy[nlocal][3];
  imageint imagecpy[nlocal];
  int atom_order[nlocal];
  int molcpy[nlocal];
  //outfile << "ch1: " << ch1 << ", ch2: " << ch2 << std::endl;
  int count=0;
  for (int i = 0; i<nlocal; i++)
  {
    //outfile << "atom_order[" << tag[i]-1 << "]: " << i << std::endl;
    atom_order[tag[i]-1] = i;
    xcpy[tag[i]-1][0] = x[i][0];
    xcpy[tag[i]-1][1] = x[i][1];
    xcpy[tag[i]-1][2] = x[i][2];
    molcpy[i] = molecule[i];
    imagecpy[tag[i]-1] = image[i];
    domain->unmap(xcpy[tag[i]-1], image[i]);
  }

  /*for (int i = 0; i<nlocal; i++)
  {
    std::cout << "atom_order[" << i << "]: " << atom_order[i] << std::endl;
    std::cout << "xcpy: " << xcpy[i][0] << ", " << xcpy[i][1] << ", " << xcpy[i][2] << std::endl;
    std::cout << "molecule[" << i << "]: " << molecule[i] << std::endl;
    std::cout << "molcpy[" << i << "]: " << molcpy[i] << std::endl;
    std::cout << "num_bond[" << i << "]: " << num_bond[i] << std::endl;
    std::cout << "bond_atom[" << atom_order[i] << "][0]: " << bond_atom[atom_order[i]][0] << std::endl;
    if (i+1 < nlocal)
    {
      std::cout << "bond_type[" << i << "][" << i+1 << "]: " << bond_type[i][i+1] << std::endl;
    }
    std::cout << std::endl;
  }*/
  //std::cout << "num bond " << num_bond[5] << std::endl;
  for (int i=0; i<nlocal-1; ++i)
  {
    if (molcpy[i] != ch1 || molcpy[i+1] != ch1) continue;
    //int i1=atom_order[i];
    //if (i+1 == nlocal || molecule[i] != molecule[i+1]) continue;
    //int i2 = (int) bond_atom[atom_order[i]][0] - 1;
    //if (i2 <= i1) continue;
    for (int j=0; j<nlocal; ++j)
    {
      if (ch1 == ch2 && j < i+1 && !pflag) continue;
      if (molcpy[j] != ch2) continue;
      if (j+1 == nlocal || molcpy[j] != molcpy[j+1]) continue;
      if (ch1 == ch2 && (i+1 == j || i == j+1 || i == j || i+1 == j+1) && (jxtrans == jytrans == jztrans == 0)) continue;
      
      lk+=compute_one(xcpy[i],xcpy[i+1],xcpy[j],xcpy[j+1], jxtrans, jytrans, jztrans);
      if (ch1 == ch2 && pflag) count++;
    }
  }

  if (ch1 == ch2 && !pflag)
  {
    lk = lk*2;
  }
  return lk;
}

double ComputeLK::compute_one(const double *x1, const double *x2, const double *x3, const double *x4, double jxtrans, double jytrans, double jztrans) {
                //method computes lk between two given edges, second edge can be translated by given amount
  a1[0]=x1[0];
  a1[1]=x1[1];
  a1[2]=x1[2];
  a2[0]=x2[0];
  a2[1]=x2[1];
  a2[2]=x2[2];
  b1[0]=x3[0]+jxtrans;
  b1[1]=x3[1]+jytrans;
  b1[2]=x3[2]+jztrans;
  b2[0]=x4[0]+jxtrans;
  b2[1]=x4[1]+jytrans;
  b2[2]=x4[2]+jztrans;

  MathExtra::sub3(a2, a1, ra);
  MathExtra::sub3(b2, b1, rb);
  MathExtra::sub3(a1, b1, r00);
  MathExtra::sub3(a1, b2, r01);
  MathExtra::sub3(a2, b1, r10);
  MathExtra::sub3(a2, b2, r11);

  MathExtra::cross3(r00, r01, v1);
  MathExtra::normalize3(v1, u1);
  MathExtra::cross3(r01, r11, v2);
  MathExtra::normalize3(v2, u2);
  MathExtra::cross3(r11, r10, v3);
  MathExtra::normalize3(v3, u3);
  MathExtra::cross3(r10, r00, v4);
  MathExtra::normalize3(v4, u4);

  d1=MathExtra::dot3(u1, u2);
  d2=MathExtra::dot3(u2, u3);
  d3=MathExtra::dot3(u3, u4);
  d4=MathExtra::dot3(u4, u1);

  as1=asin(d1);
  as2=asin(d2);
  as3=asin(d3);
  as4=asin(d4);

  MathExtra::cross3(ra, rb, aux1);
  aux=MathExtra::dot3(aux1, r00);
  if (aux<0) {alk=-1*((as1+as2+as3+as4)/(4*MathConst::MY_PI));}
  else {alk=(as1+as2+as3+as4)/(4*MathConst::MY_PI);}
  if (isnan(alk)) alk=0;
  return alk;
}