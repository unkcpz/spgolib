#ifndef _WRAPPER_H
#define _WRAPPER_H
#endif

#include "spglib.h"

int spgo_delaunay_reduce(double lattice[], const double symprec);
int spgo_standardize_cell(double lattice[],
                        double position[],
                        int types[],
                        const int num_atom,
                        const int to_primitive,
                        const int no_idealize,
                        const double symprec);
SpglibDataset * spgo_get_dataset(double lattice[],
                                 double position[],
                                 int types[],
                                 const int num_atom,
                                 const double symprec);
void spgo_free_dataset(SpglibDataset *dataset);
int spgo_spacegroup_number(SpglibDataset *dataset);
int spgo_hall_number(SpglibDataset *dataset);
int spgo_international_symbol(SpglibDataset *dataset, char out[11]);
int spgo_hall_symbol(SpglibDataset *dataset, char out[17]);
int spgo_dataset_n_operations(SpglibDataset *dataset);
int spgo_dataset_rotations(SpglibDataset *dataset, int *out);
int spgo_dataset_tranlations(SpglibDataset *dataset, double *out);

void flat_mat_3D(double mat[][3], double flat[], int n);
void mat_flat_3D(double flat[], double mat[][3], int n);
