#include "spglib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static void test_spg_find_primitive_BCC(void);
static void test_spg_standardize_cell_BCC_prim(void);
static int sub_spg_standardize_cell(double lattice[3][3],
				    double position[][3],
				    int types[],
				    const int num_atom,
				    const double symprec,
				    const int to_primitive,
				    const int no_idealize);
int spgo_find_primitive(double lattice[],
                        double position[],
                        int types[],
                        const int num_atom,
                        const double symprec);
static void show_cell(double lattice[3][3],
		      double position[][3],
		      const int types[],
		      const int num_atom);
static void showgo_cell(double lattice[],
		      double position[],
		      const int types[],
		      const int num_atom);

int main(void)
{
  test_spg_find_primitive_BCC();
  // test_spg_standardize_cell_BCC_prim();

  return 0;
}

static void test_spg_find_primitive_BCC(void)
{
  double lattice[] = {4, 0, 0, 0, 4, 0, 0, 0, 4};
  double position[] = {
    0, 0, 0,
    0.5, 0.0, 0.5,
    0.0, 0.5, 0.5,
    0.5, 0.5, 0.0,
  };
  int types[] = {1, 1, 1, 1};
  int i, num_atom = 4, num_primitive_atom;
  double symprec = 1e-5;

  /* lattice, position, and types are overwirtten. */
  printf("*** Example of spg_find_primitive (BCC unitcell --> primitive) ***:\n");
  showgo_cell(lattice, position, types, num_atom);
	printf("num_atom: %d, num_primitive_atom: %d\n", num_atom, num_primitive_atom);
  num_primitive_atom = spgo_find_primitive(lattice, position, types, num_atom, symprec);
  if (num_primitive_atom == 0) {
    printf("Primitive cell was not found.\n");
  } else {
    showgo_cell(lattice, position, types, num_primitive_atom);
  }
}

int spgo_find_primitive(double lattice[],
                        double position[],
                        int types[],
                        const int num_atom,
                        const double symprec) {

	double latt[3][3];
	double pos[num_atom][3];
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			latt[i][j] = lattice[3*i+j];
		}
	}
	for(int i=0; i<num_atom; i++){
		for(int j=0; j<3; j++){
			pos[i][j] = position[3*i+j];
		}
	}
	int r=0;
  r = spg_find_primitive(latt, pos, types, num_atom, symprec);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			lattice[3*i+j] = latt[i][j];
		}
	}
	for(int i=0; i<num_atom; i++){
		for(int j=0; j<3; j++){
			position[3*i+j] = pos[i][j];
		}
	}

	return r;
}

static void test_spg_standardize_cell_BCC_prim(void)
{
  double lattice[3][3] = {{-2.01, 2, 2}, {2, -2.02, 2}, {2, 2, -2.03}};
  double position[][3] = {
    {0.002, 0, 0},
  };
  int types[] = {1};
  int i, j, k, num_atom = 1, num_primitive_atom;
  double symprec = 1e-1;

  /* lattice, position, and types are overwirtten. */
  printf("*** Example of spg_standardize_cell (BCC primitive) ***:\n");
  printf("------------------------------------------------------\n");
  for (j = 0; j < 2; j++) {
    for (k = 0; k < 2; k++) {
      sub_spg_standardize_cell(lattice,
			       position,
			       types,
			       num_atom,
			       symprec,
			       j,
			       k);
      printf("------------------------------------------------------\n");
    }
  }
}

static int sub_spg_standardize_cell(double lattice[3][3],
				    double position[][3],
				    int types[],
				    const int num_atom,
				    const double symprec,
				    const int to_primitive,
				    const int no_idealize)
{
  int i, num_primitive_atom;
  double lat[3][3], pos[num_atom][3];
  int typ[num_atom];

  for (i = 0; i < 3; i++) {
    lat[i][0] = lattice[i][0];
    lat[i][1] = lattice[i][1];
    lat[i][2] = lattice[i][2];
  }

  for (i = 0; i < num_atom; i++) {
    pos[i][0] = position[i][0];
    pos[i][1] = position[i][1];
    pos[i][2] = position[i][2];
    typ[i] = types[i];
  }

  /* lattice, position, and types are overwirtten. */
  num_primitive_atom = spg_standardize_cell(lat,
					    pos,
					    typ,
					    num_atom,
					    to_primitive,
					    no_idealize,
					    symprec);
  printf("VASP POSCAR format: ");
  if (to_primitive == 0) {
    printf("to_primitive=0 and ");
  } else {
    printf("to_primitive=1 and ");
  }

  if (no_idealize == 0) {
    printf("no_idealize=0\n");
  } else {
    printf("no_idealize=1\n");
  }
  printf("1.0\n");
  for (i = 0; i < 3; i++) {
    printf("%f %f %f\n", lat[0][i], lat[1][i], lat[2][i]);
  }
  printf("%d\n", num_primitive_atom);
  printf("Direct\n");
  for (i = 0; i < num_primitive_atom; i++) {
    printf("%f %f %f\n", pos[i][0], pos[i][1], pos[i][2]);
  }
}
static void show_cell(double lattice[3][3],
		      double position[][3],
		      const int types[],
		      const int num_atom)
{
  int i;

  printf("Lattice parameter:\n");
  for (i = 0; i < 3; i++) {
    printf("%f %f %f\n", lattice[0][i], lattice[1][i], lattice[2][i]);
  }
  printf("Atomic positions:\n");
  for (i = 0; i < num_atom; i++) {
    printf("%d: %f %f %f\n",
	   types[i], position[i][0], position[i][1], position[i][2]);
  }
}

static void showgo_cell(double lattice[],
		      double position[],
		      const int types[],
		      const int num_atom)
{
  int i;
	double latt[3][3];
	double pos[num_atom][3];
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			latt[i][j] = lattice[3*i+j];
		}
	}
	for(int i=0; i<num_atom; i++){
		for(int j=0; j<3; j++){
			pos[i][j] = position[3*i+j];
		}
	}

  printf("Lattice parameter:\n");
  for (i = 0; i < 3; i++) {
    printf("%f %f %f\n", latt[0][i], latt[1][i], latt[2][i]);
  }
  printf("Atomic positions:\n");
  for (i = 0; i < num_atom; i++) {
    printf("%d: %f %f %f\n",
	   types[i], pos[i][0], pos[i][1], pos[i][2]);
  }
}
