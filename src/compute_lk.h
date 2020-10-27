#ifdef COMPUTE_CLASS

ComputeStyle(lk, ComputeLK)

#else

#ifndef LMP_COMPUTE_LK_H
#define LMP_COMPUTE_LK_H

#include "compute.h"

namespace LAMMPS_NS {

  class ComputeLK : public Compute {
  public:
    ComputeLK(class LAMMPS *, int, char **);
    ~ComputeLK();
    void init();
    void compute_array();
    void compute_vector();
    void init_list(int, class NeighList *);
    

  private:
    void allocate();
    double compute_one(const double *, const double *, const double *, const double *, double, double, double);
    double lk(int, int, double=0, double=0, double=0);
    class FixStore *fix;
    class NeighList *list;
    tagint *frmolecule;
    bool pflag;
    bool grflag;
    bool frflag;
    bool wrflag;
    int jgroup,jgroupbit;
    char *group2;
    char *id_fix;
    int ix, iy, iz, jx, jy, jz;
    int nchains;
    int nmols;
    int *mols;
    bool found;
    int current;
    int atom1;
    int atom2;
    int atom3;
    int atom4;
    int nvec;
    int offset;
    int lastid;
    int nstates;
    int currentid;
    double torsion;
    double scale;
    int point1;
    int point2;
    int point3;
    int point4;
    int o1;
    int o2;
    int o3;
    int o4;
    bool os1;
    bool os2;
    bool os3;
    bool os4;
    bool do_linking;
    char *idfragment;
    char *idproperty;
    bool samechain;
    int m;
    int first;
    class ComputeFragmentAtom *cfragment;
    class ComputePropertyAtom *cproperty;
    double a1 [3];
    double a2 [3];
    double b1 [3];
    double b2 [3];
    double c1 [3];
    double sum;
    double ra [3];
    double rb [3];
    double r00 [3];
    double r01 [3];
    double r10 [3];
    double r11 [3];
    double v1 [3];
    double v2 [3];
    double v3 [3];
    double v4 [3];
    double u1 [3];
    double u2 [3];
    double u3 [3];
    double u4 [3];
    double d1;
    double d2;
    double d3;
    double d4;
    double as1;
    double as2;
    double as3;
    double as4;
    double aux1 [3];
    double aux;
    double alk;
    double norm;
    double sign;
    tagint mol1ID;
    tagint mol2ID;
    int nmol1;
    int nmol2;
    int count1;
    int count2;
    double acn;
    double maximum;
    int maxch1;
    int maxch2;
    double maxxtrans;
    double maxytrans;
    double maxztrans;
  };

}

#endif
#endif
