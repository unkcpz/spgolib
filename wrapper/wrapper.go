package wrapper

//#cgo CFLAGS: -g
//#cgo LDFLAGS: -lm
// #include <stdlib.h>
//#include "wrapper.h"
import "C"
import "unsafe"

type SpglibDataset struct {
  // Spacegrou number and symbol
  SpacegroupNumber int
  SpacegroupSymbol string
  // Hall number and symbol
  HallNumber int
  HallSymbol string
  // Pointgroup symbol
  PointgroupSymbol string

  // wyckoffs
  Wyckoffs []int
  // equivalent atoms
  AtomEquivalent []int

  // Noperations number of operations
  Noperations int
  // Rotations are 1D represent of rotations tobe reshape to n*(3*3)
  Rotations []int
  // Translations are 1D represent of tarnslations tobe reshape to n*(3)
  Translations []float64
}

// NewSpglibDataset create SpglibDataset with a cell
func NewSpglibDataset(
  lattice []float64,
  position []float64,
  types []int,
  num_atom int,
  symprec float64,
  ) *SpglibDataset {

  c_dataset := getDataset(lattice, position, types, num_atom, symprec)
  defer freeDataset(c_dataset)

  spacegroupNumber := C.spgo_spacegroup_number(c_dataset)

  ptr := C.malloc(C.sizeof_char * 11)
  size := C.spgo_international_symbol(c_dataset, (*C.char)(ptr))
  b := C.GoBytes(ptr, size)
  C.free(unsafe.Pointer(ptr))
  internationalSymbol := string(b)

  hallNumber := C.spgo_hall_number(c_dataset)

  ptr = C.malloc(C.sizeof_char * 17)
  size = C.spgo_hall_symbol(c_dataset, (*C.char)(ptr))
  b = C.GoBytes(ptr, size)
  C.free(unsafe.Pointer(ptr))
  hallSymbol := string(b)

  ptr = C.malloc(C.sizeof_char * 6)
  size = C.spgo_pointgroup_symbol(c_dataset, (*C.char)(ptr))
  b = C.GoBytes(ptr, size)
  C.free(unsafe.Pointer(ptr))
  pointgroupSymbol := string(b)

  nOp := C.spgo_dataset_n_operations(c_dataset)

  rotations := make([]int32, nOp*9, nOp*9)
  C.spgo_dataset_rotations(c_dataset, (*C.int)(&rotations[0]))

  translations := make([]float64, nOp*3, nOp*3)
  C.spgo_dataset_tranlations(c_dataset, (*C.double)(&translations[0]))

  wyckoffs := make([]int32, num_atom)
  C.spgo_dataset_wyckoffs(c_dataset, (*C.int)(&wyckoffs[0]))

  equivalent_atoms := make([]int32, num_atom)
  C.spgo_dataset_equivalent_atoms(c_dataset, (*C.int)(&equivalent_atoms[0]))

  return &SpglibDataset{
    SpacegroupNumber: int(spacegroupNumber),
    SpacegroupSymbol:  internationalSymbol,
    HallNumber: int(hallNumber),
    HallSymbol: hallSymbol,
    PointgroupSymbol: pointgroupSymbol,
    Wyckoffs: tobit(wyckoffs),
    AtomEquivalent: tobit(equivalent_atoms),
    Noperations: int(nOp),
    Rotations: tobit(rotations),
    Translations: translations,
  }
}

func getDataset(
  lattice []float64,
  position []float64,
  types []int,
  num_atom int,
  symprec float64) *C.SpglibDataset {

  tp32 := to32bit(types)
  return C.spgo_get_dataset(
    (*C.double)(unsafe.Pointer(&lattice[0])),
    (*C.double)(unsafe.Pointer(&position[0])),
    (*C.int)(unsafe.Pointer(&tp32[0])),
    (C.int)(num_atom),
    (C.double)(symprec),
  )
}

func freeDataset(dataset *C.SpglibDataset) {
  C.spgo_free_dataset(dataset)
}

// DelaunayReduce reduce the cell by delaunay approach
func DelaunayReduce(lattice []float64, symprec float64) []float64 {
  ret := C.spgo_delaunay_reduce(
    (*C.double)(unsafe.Pointer(&lattice[0])),
    (C.double)(symprec),
  )
  if ret == 0 {
    panic("spg_delaunay_reduce failed")
  }
  return lattice
}

// StandardizeCell standardize cell 
func StandardizeCell(
  lattice []float64,
  position []float64,
  types []int,
  num_atom int,
  to_primitive int,
  no_idealize int,
  symprec float64) ([]float64, []float64, []int) {

  tp32 := to32bit(types)

  originL := make([]float64, len(lattice))
  copy(originL, lattice)
  originP := make([]float64, len(position))
  copy(originP, position)
  originT := make([]int32, len(tp32))
  copy(originT, tp32)

  n := C.spgo_standardize_cell(
    (*C.double)(unsafe.Pointer(&lattice[0])),
    (*C.double)(unsafe.Pointer(&position[0])),
    (*C.int)(unsafe.Pointer(&tp32[0])),
    (C.int)(num_atom),
    (C.int)(to_primitive),
    (C.int)(no_idealize),
    (C.double)(symprec),
  )
  if n == 0 {
    panic("spg_standardize_cell failed")
  }
  // Recalculate slice dynamiclly
  if int(n) > cap(position)/3 {
    outP := make([]float64, n*3, n*3)
    for i:=0; i<len(originP); i++ {
      outP[i] = originP[i]
    }
    outT := make([]int32, n, n)
    for i:=0; i<len(originT); i++ {
      outT[i] = originT[i]
    }
    C.spgo_standardize_cell(
      (*C.double)(unsafe.Pointer(&originL[0])),
      (*C.double)(unsafe.Pointer(&outP[0])),
      (*C.int)(unsafe.Pointer(&outT[0])),
      (C.int)(num_atom),
      (C.int)(to_primitive),
      (C.int)(no_idealize),
      (C.double)(symprec),
    )
    tp := tobit(outT)
    return originL, outP, tp
  }

  tp := tobit(tp32)
  return lattice, position, tp
}

// C.int is 32 bit even on a 64bit system, but Go int is 32 or 64 bit.
// So we need to convert in order to pass C int arrays.
func to32bit(a []int) []int32 {
	b := make([]int32, len(a))
	for i := range b {
		b[i] = int32(a[i])
	}
	return b
}

func tobit(a []int32) []int {
	b := make([]int, len(a))
	for i := range b {
		b[i] = int(a[i])
	}
	return b
}
