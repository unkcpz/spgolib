package wrapper

import (
  "testing"
  // "fmt"
)

func TestNewSpglibDatabase(t *testing.T) {
  lattice := []float64{
    4, 0, 0,
    0, 4, 0,
    0, 0, 4,
  }
  position := []float64{
    0, 0, 0,
    0.5, 0.5, 0.5,
  }
  types := []int{1, 1}
  n_atoms := 2

  goSpglibDataset := NewSpglibDataset(
    lattice,
    position,
    types,
    n_atoms,
    1e-5,
  )

  got_0 := goSpglibDataset.SpacegroupNumber
  if got_0 != 229 {
    t.Errorf("spacegroupNumber, expect 229, got %d", got_0)
  }

  got_1 := goSpglibDataset.InternationalSymbol
  if got_1 != "Im-3m" {
    t.Errorf("international_symbol, expect Im-3m, got %s", got_1)
  }

  got_2 := goSpglibDataset.HallNumber
  if got_2 != 529 {
    t.Errorf("hallNumber, expect 229, got %d", got_2)
  }

  got_3 := goSpglibDataset.HallSymbol
  if got_3 != "-I 4 2 3" {
    t.Errorf("international_symbol, expect '-I 4 2 3', got %s", got_3)
  }

  got_4 := goSpglibDataset.Noperations
  if got_4 != 96 {
    t.Errorf("number of operations, expect 96, got %d", got_4)
  }

  // fmt.Println(goSpglibDataset.Rotations)
  // fmt.Println(goSpglibDataset.Translations)
}

func TestDelaunayReduce(t *testing.T) {
  lattice := []float64{4, 0, 0, 2, 4, 0, 0, 0, 4}
  origin := make([]float64, len(lattice))
  copy(origin, lattice)

  lattice = DelaunayReduce(lattice, 1e-5)

  expected := []float64{0, 0, -4, 0, -4, -2, -4, 0, 0}
  for i, _ := range lattice {
    if lattice[i] - expected[i] > 1e-5 {
      t.Errorf("delaunay reduce lattice: %v, expected %v, got %v.", origin, expected, lattice)
      break
    }
  }
}

func TestStandardizeCellBCCRefine(t *testing.T) {
  // BCC to its Conventional
  lattice := []float64{
    -2, 2, 2,
    2, -2, 2,
    2, 2, -2,
  }
  originL := make([]float64, len(lattice))
  copy(originL, lattice)

  position := []float64{
    0, 0, 0,
  }
  originP := make([]float64, len(position))
  copy(originP, position)

  types := []int{1}
  n_atoms := 1
  lattice, position, types = StandardizeCell(lattice, position, types, n_atoms, 0, 0, 1e-5)

  expectedL := []float64{4, 0, 0, 0, 4, 0, 0, 0, 4}
  expectedP := []float64{0, 0, 0, 0.5, 0.5, 0.5}
  expectedT := []int{1, 1}
  for i, _ := range expectedL {
    if lattice[i] - expectedL[i] > 1e-5 {
      t.Error("StandardizeCell error")
      break
    }
  }
  for i, _ := range expectedP{
    if position[i] - expectedP[i] > 1e-5 {
      t.Error("StandardizeCell error")
      break
    }
  }
  for i, _ := range expectedT {
    if types[i] != expectedT[i] {
      t.Error("StandardizeCell error")
      break
    }
  }
}

func TestStandardizeCellBCCPrimitive(t *testing.T) {
  // BCC to its primitive
  lattice := []float64{
    4, 0, 0,
    0, 4, 0,
    0, 0, 4,
  }
  originL := make([]float64, len(lattice))
  copy(originL, lattice)

  position := []float64{
    0, 0, 0,
    0.5, 0.5, 0.5,
  }
  originP := make([]float64, len(position))
  copy(originP, position)

  types := []int{1, 1}
  n_atoms := 2
  lattice, position, types = StandardizeCell(lattice, position, types, n_atoms, 1, 1, 1e-5)

  expectedL := []float64{-2, 2, 2, 2, -2, 2, 2, 2, -2}
  expectedP := []float64{0, 0, 0}
  expectedT := []int{1}
  for i, _ := range expectedL {
    if lattice[i] - expectedL[i] > 1e-5 {
      t.Errorf("got %v, expected %v", lattice, expectedL)
      break
    }
  }
  for i, _ := range expectedP{
    if position[i] - expectedP[i] > 1e-5 {
      t.Error("StandardizeCell error")
      break
    }
  }
  for i, _ := range expectedT {
    if types[i] != expectedT[i] {
      t.Error("StandardizeCell error")
      break
    }
  }
}
