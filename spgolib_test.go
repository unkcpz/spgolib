package spgolib

import (
	// "fmt"
	"testing"
)

func TestNewDataset(t *testing.T) {
	ds := NewDataset(
		[]float64{2, 0, 0, 1, 1.732050, 0, 0, 0, 14},
		[]float64{0, 0, 0},
		[]int{1},
		1e-5)

	got_int := ds.SpaceNumber
	if got_int != 191 {
		t.Errorf("spacegroupNumber, expect 229, got %d", got_int)
	}

	got_str := ds.SpaceSymbol
	if got_str != "P6/mmm" {
		t.Errorf("international_symbol, expect Im-3m, got %s", got_str)
	}

	got_int = ds.HallNumber
	if got_int != 485 {
		t.Errorf("hallNumber, expect 485, got %d", got_int)
	}

	got_str = ds.HallSymbol
	if got_str != "-P 6 2" {
		t.Errorf("international_symbol, expect '-P 6 2', got %s", got_str)
	}

	got_str = ds.PointSymbol
	if got_str != "6/mmm" {
		t.Errorf("international_symbol, expect '6/mmm', got %s", got_str)
	}

	got_int = ds.Nops
	if got_int != 24 {
		t.Errorf("number of operations, expect 24, got %d", got_int)
	}
	// fmt.Println(ds.Rotations)
	// fmt.Println(ds.Translations)
}

func TestDelaunayReduce(t *testing.T) {
	lattice := []float64{4, 0, 0, 2, 4, 0, 0, 0, 4}
	origin := make([]float64, len(lattice))
	copy(origin, lattice)

	lattice = Delaunay(lattice, 1e-5)

	expected := []float64{-4, 0, 2, 0, 0, -4, 0, -4, 0}
	for i, _ := range lattice {
		if lattice[i]-expected[i] > 1e-5 {
			t.Errorf("delaunay reduce lattice: %v, expected %v, got %v.", origin, expected, lattice)
			break
		}
	}
}

func TestStandardizeBCCRefine(t *testing.T) {
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
	lattice, position, types = Standardize(lattice, position, types, n_atoms, false, false, 1e-5)

	expectedL := []float64{4, 0, 0, 0, 4, 0, 0, 0, 4}
	expectedP := []float64{0, 0, 0, 0.5, 0.5, 0.5}
	expectedT := []int{1, 1}
	for i, _ := range expectedL {
		if lattice[i]-expectedL[i] > 1e-5 {
			t.Error("StandardizeCell error")
			break
		}
	}
	for i, _ := range expectedP {
		if position[i]-expectedP[i] > 1e-5 {
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

func TestStandardizeBCCPrimitive(t *testing.T) {
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
	lattice, position, types = Standardize(lattice, position, types, n_atoms, true, true, 1e-5)

	expectedL := []float64{-2, 2, 2, 2, -2, 2, 2, 2, -2}
	expectedP := []float64{0, 0, 0}
	expectedT := []int{1}
	for i, _ := range expectedL {
		if lattice[i]-expectedL[i] > 1e-5 {
			t.Errorf("got %v, expected %v", lattice, expectedL)
			break
		}
	}
	for i, _ := range expectedP {
		if position[i]-expectedP[i] > 1e-5 {
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
