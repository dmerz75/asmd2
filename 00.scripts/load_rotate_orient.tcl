# #!/usr/bin/tclsh

# load_rotate_orient.tcl
# This script orients a protein on the Z-axis.

# mol new
# mol load, mol list, mol delete 7

# python equivalents
# load('pdb','alanine.pdb')
# load('psf','alanine.psf')

# set molecule [lindex $argv 0]
# puts $molecule
# exit

# source from TKCON
proc orient_protein_on_z { molec start end } {
    # protein
    set all [atomselect $molec all]

    set resid_CA_begin_coord [atomselect $molec "resid $start and name CA"]
    set resid_CA_end_coord [atomselect $molec "resid $end and name CA"]

    set pos_begin [measure center $resid_CA_begin_coord]
    set inv_pos_begin [vecinvert [measure center $resid_CA_begin_coord]]
    set pos_end   [measure center $resid_CA_end_coord]

    # get the protein vector; $end - $begin for vecsub
    set protein_vector [vecsub $pos_end $pos_begin]
    set prot_on_x [transvecinv $protein_vector]

    # translate from x0,y0,z0 to {0.1 -0.5 0.08}
    set inv_protein_vector [vecinvert [vecsub $pos_end $pos_begin]]
    set inv_protein_vector [vecadd $inv_protein_vector {0.1 -0.05 0.08}]

    # translate: (1) near origin (2) to x axis
    # (1)
    # $all moveby $inv_pos_begin
    # (2)
    $all move $prot_on_x
    $all move [trans y -90]
    $all moveby {0.1 -0.5 0.08}
}

proc draw_origin { molec x y z} {
    # draw dashed red lines on the scaled unit vectors
    puts "$molec"
    # set mol_id [mol new]
    # puts "$scaled_vector"
    set scaled_vector [vecscale $z {0 0 1}]
    graphics $molec color blue
    graphics $molec line {0 0 0} $scaled_vector width 8 style dashed
    set scaled_vector [vecscale $y {0 1 0}]
    graphics $molec color green
    graphics $molec line {0 0 0} $scaled_vector width 8 style dashed
    set scaled_vector [vecscale $x {1 0 0}]
    graphics $molec color red
    graphics $molec line {0 0 0} $scaled_vector width 8 style dashed
    # mol delete $molec
}

# Executed Commands
# argv molec first_residue second_residue
orient_protein_on_z 0 1 10
# argv molec scaling
set mol_id [mol new]
draw_origin $mol_id 5.0 5.0 15.0

# set outfile [open rmsd.dat w];                                             
# set nf [molinfo top get numframes]
# set frame0 [atomselect top "protein and backbone and noh" frame 0]
# set sel [atomselect top "protein and backbone and noh"]
# # rmsd calculation loop
# for {set i 1 } {$i < $nf } { incr i } {
#     $sel frame $i
#     $sel move [measure fit $sel $frame0]
#     puts $outfile "[measure rmsd $sel $frame0]"
# }
# close $outfile
