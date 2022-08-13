
if {1} {
	foreach X { A B D E H L P Y} {
		if { [expr {![catch {file lstat PMEPUFA/${X} finfo}]}] == 0 } {	
			exec cp -r PMEPUFA/XXX PMEPUFA/${X}	
		}
	}
}

set L {}

foreach i {A B D E H L P Y} {
	set $i [mol new PMEPUFA/equb3/dcdb3/${i}.psf]
	mol addfile PMEPUFA/equb3/dcdb3/${i}.1.dcd waitfor all
	
	set sel [atomselect [set $i] "not water and not resname CLA POT"]
	set m [measure center $sel]
	foreach {a b c} $m {break}
	set d [expr -1.0*$c]
	set e [subst {0 0 $d}]	
	$sel moveby $e
	
	$sel writepdb temp.pdb
	$sel writepsf temp.psf
	
	resetpsf
	readpsf temp.psf
	readpsf PME.psf
	coordpdb temp.pdb
	coordpdb PME.pdb	

	writepsf tempr.psf
	writepdb tempr.pdb

	solvate tempr.psf tempr.pdb -minmax {{-32.5 -32.5 -58} {32.5 32.5 -21}} -b 2 -o temps
	solvate temps.psf temps.pdb -minmax {{-32.5 -32.5 21} {32.5 32.5 58}} -s WA -b 2 -o tempt
	
	autoionize -psf tempt.psf -pdb  tempt.pdb -sc 0.15 -cation POT -o tempu

	mol delete all
	mol new tempu.psf
	mol addfile tempu.pdb
	


	###
	### BOUNDARY PHOSPHO
	###
	
	set all [atomselect top all]
	$all set beta 0
	$all set occupancy 0
	$all writepsf PMEPUFA/${i}/${i}.psf	
	$all writepdb PMEPUFA/${i}/${i}.pdb

	set sel [atomselect top "name P P1 P2 and not ( sqr(x)+sqr(y) < sqr(25) )"]	
	$sel set beta 1
	$all writepdb PMEPUFA/${i}/${i}_border.ref 	
	$all set beta 0	
	
	

	###
	### TAIL 
	###

	set sel [atomselect top "resname 6FH and name CAA"]
	$sel set occupancy 1 
	$all writepdb PMEPUFA/${i}/${i}_PME.ref
	$all set occupancy 0	
	$sel set beta 1
	$all writepdb PMEPUFA/${i}/${i}_centerfrom.ref

	###
	### CENTER
	###

	$all set beta 0
	set sel [atomselect top "sqr(x) + sqr(y) + sqr(z - 14) < sqr(5)"]
	$sel set beta 1
	$all writepdb PMEPUFA/${i}/${i}_centerto.ref
	
}

###modified after moving to equ
if {1} {
	foreach X {A B D E H L P Y} {
		exec sed -i "s/XXX/${X}/g" PMEPUFA/${X}/colvar.str
		exec sed -i "s/XXX/${X}/g" PMEPUFA/${X}/s0.inp
		exec sed -i "s/XXX/${X}/g" PMEPUFA/${X}/s0.qsub

 		for {set i 1} {$i <= 3} {incr i} {
			exec sed "s/restart 0/restart ${i}/g" PMEPUFA/${X}/s0.inp > PMEPUFA/${X}/s${i}.inp
			exec sed "s/s0/s${i}/g" PMEPUFA/${X}/s0.qsub > PMEPUFA/${X}/s${i}.qsub
			exec sed -i "s/${X}.0/${X}.${i}/g" PMEPUFA/${X}/s${i}.qsub
		}			
	}
}

mol delete all 

if {1} {
	foreach i {A B D E H L P Y} {
		mol new PMEPUFA/${i}/${i}.psf
		mol addfile PMEPUFA/${i}/${i}.pdb
		mol off top
	}
}


if {1} {
	foreach i {A B D E H L P Y} {
		mol new PMEPUFA/${i}/${i}_border.ref type pdb
		mol modselect 0 top "beta 1"
		mol modstyle 0 top CPK	
		set sel [atomselect top "beta 1"]
		puts [$sel num]	
	}
}

if {1} {
	foreach i {A B D E H L P Y} {
		mol new PMEPUFA/${i}/${i}_PME.ref type pdb
		mol modselect 0 top "occupancy 1"
		mol modstyle 0 top CPK		
	}
}

if {1} {
	foreach i [glob *temp*] { exec rm $i }
	#foreach i [glob *tmp*] { exec rm $i }
}

puts $L
