package spgolib

import (
	"github.com/unkcpz/spgolib/wrapper"
	"gonum.org/v1/gonum/mat"
)

type Dataset struct {
	// Spacegrou number and symbol
	SpaceNumber int
	SpaceSymbol string
	// Hall number and symbol
	HallNumber int
	HallSymbol string
	// Pointgroup symbol
	PointSymbol string

	// wyckoffs
	Wyckoffs []int
	// equivalent atoms
	AtomEquivalent []int

	// Noperations number of operations
	Nops int
	// Rotations are reshaped rotations n*(3*3)
	Rotations []*mat.Dense
	// Translations are reshaped tarnslations n*(3)
	Translations []*mat.VecDense
}

// NewDataset create spglib dataset from lattice (row present) position and types
func NewDataset(
	lattice,
	positions []float64,
	types []int,
	eps float64) *Dataset {

	colLattice := transpose(lattice)
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
	for i := 0; i < nops; i++ {
		rots[i] = mat.NewDense(3, 3, intToFloat(ds.Rotations[i*9:i*9+9]))
		trans[i] = mat.NewVecDense(3, ds.Translations[i*3:i*3+3])
	}
	rds := &Dataset{
		SpaceNumber:    ds.SpacegroupNumber,
		SpaceSymbol:    ds.SpacegroupSymbol,
		HallNumber:     ds.HallNumber,
		HallSymbol:     ds.HallSymbol,
		PointSymbol:    ds.PointgroupSymbol,
		Wyckoffs:       ds.Wyckoffs,
		AtomEquivalent: ds.AtomEquivalent,
		Nops:           ds.Noperations,
		Rotations:      rots,
		Translations:   trans,
	}
	return rds
}

// Delaunay reduce lattice by delaunay approach
func Delaunay(lattice []float64, symprec float64) []float64 {
	colLattice := transpose(lattice)
	return wrapper.DelaunayReduce(colLattice, symprec)
}

// Standardize return the primitive(toPrimitive=true) or conventional(toPrimitive=false) cell
func Standardize(
	lattice,
	positions []float64,
	types []int,
	Natoms int,
	toPrimitive,
	noIdealize bool,
	symprec float64) ([]float64, []float64, []int) {

	colLattice := transpose(lattice)

	var to_primitive int
	var no_idealize int
	if toPrimitive {
		to_primitive = 1
	}
	if noIdealize {
		no_idealize = 1
	}

	return wrapper.StandardizeCell(
		colLattice,
		positions,
		types,
		Natoms,
		to_primitive,
		no_idealize,
		symprec,
	)
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
