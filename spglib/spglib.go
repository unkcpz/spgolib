package spglib

import (

  "gonum.org/v1/gonum/mat"
  "github.com/unkcpz/spglib/go/wrapper"
)

type Dataset struct {
  SpaceNumber int
  SpaceSymbol string
  HallNumber int
  HallSymbol string

  Nops int
  Rotations []*mat.Dense
  Translations []*mat.VecDense
}

func NewDataset(
  rowLattice []float64,
  positions []float64,
  types []int,
  eps float64) *Dataset {

  colLattice := transpose(rowLattice)
  ds := wrapper.NewSpglibDataset(
    colLattice,
    positions,
    types,
    len(types),
    eps,
  )

  nops := ds.Noperations
  rots := make([]*mat.Dense, nops, nops)
  trans := make([]*mat.VecDense, nops, nops)
  for i:=0; i<nops; i++ {
    rots[i] = mat.NewDense(3, 3, intToFloat(ds.Rotations[i*9:i*9+9]))
    trans[i] = mat.NewVecDense(3, ds.Translations[i*3:i*3+3])
  }
  rds := &Dataset{
    SpaceNumber: ds.SpacegroupNumber,
    SpaceSymbol: ds.InternationalSymbol,
    HallNumber: ds.HallNumber,
    HallSymbol: ds.HallSymbol,
    Nops: ds.Noperations,
    Rotations: rots,
    Translations: trans,
  }
  return rds
}

func transpose(latt []float64) []float64 {
  r := make([]float64, 9, 9)
  rx := 0
  for _, e := range latt {
    r[rx] = e
    rx += 3
    if rx >= 9 {
      rx -= 8
    }
  }
  return r
}

func intToFloat(a []int) []float64 {
  b := make([]float64, len(a))
  for i := range b {
    b[i] = float64(a[i])
  }
  return b
}

func MatData(mat mat.Dense) []float64 {
  blasM := mat.RawMatrix()
  return blasM.Data
}
