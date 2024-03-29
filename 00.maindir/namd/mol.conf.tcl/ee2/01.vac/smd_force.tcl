proc calculateDistance {c1 c2} {
 set x1 [lindex $c1 0]
 set y1 [lindex $c1 1] 
 set z1 [lindex $c1 2]
 set x2 [lindex $c2 0] 
 set y2 [lindex $c2 1]
 set z2 [lindex $c2 2]
 return [expr sqrt(pow($x1-$x2,2) + pow($y1-$y2,2) + pow($z1-$z2,2))]
}

# Atoms selected for force application 

set id1 [atomid PEP 1 CA]
set grp1 {}
lappend grp1 $id1
set a1 [addgroup $grp1]

set id2 [atomid PEP 10 CA]
set grp2 {}
lappend grp2 $id2
set a2 [addgroup $grp2]

# set the output frequency, initialize the time counter
set Tclfreq xxfreqxx
set t 0
#set currentStep yyyyy

# contraint points

set c1x 0.0
set c1y 0.0
set c1z 0.0

set c2x 0.0
set c2y 0.0
#set c2z [expr xxzcoordxx+$currentStep*2]
#set c2z xxcur_zxx
set c2z [expr xxzcoordxx+xxcur_zxx]

# force constant (kcal/mol/A^2)
set k 7.2

# pulling velocity (A/timestep)
set v xxvelocityxx

set outfilename smdforces.out
open $outfilename w

proc calcforces {} {

  global Tclfreq t k v a1 a2 c1x c1y c1z c2x c2y c2z instantD outfilename

  # get coordinates

  loadcoords coordinate

  set r1 $coordinate($a1)
  set r1x [lindex $r1 0]
  set r1y [lindex $r1 1]
  set r1z [lindex $r1 2]

  set r2 $coordinate($a2)
  set r2x [lindex $r2 0]
  set r2y [lindex $r2 1]
  set r2z [lindex $r2 2]

  # calculate forces

  set f1x [expr $k*($c1x-$r1x)]
  set f1y [expr $k*($c1y-$r1y)]
  set f1z [expr $k*($c1z-$r1z)]
  lappend f1 $f1x $f1y $f1z

  set f2x [expr $k*($c2x-$r2x)]
  set f2y [expr $k*($c2y-$r2y)]
  set f2z [expr $k*($c2z+$v*$t-$r2z)]
  lappend f2 $f2x $f2y $f2z

  # apply forces

  addforce $a1 $f1
  addforce $a2 $f2

  # output

  set foo [expr $t % $Tclfreq]
  if { $foo == 0 } {
      set instantD [calculateDistance $r1 $r2]
      set outfile [open $outfilename a]
      set time [expr $t*xxtsxx/1000.0]
      puts $outfile "$time $r2z $instantD $f2z"
      close $outfile
  }
  incr t
  return
}
