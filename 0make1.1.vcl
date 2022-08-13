source ../0rotate.pgn
package require psfgen
resetpsf

if {1} {
	foreach X { A B D E H L P Y} {
		if { [expr {![catch {file lstat PMEPUFA/${X}b finfo}]}] == 0 } {	
			exec cp -r PMEPUFA/YYY PMEPUFA/${X}b	
		}
	}
}

set LLPOPEnum 55
set LLLPEAnum 0
set LLLPEBnum 1
set LLPOPGnum 4
set LLCLnum 2
set LLPXPEnum 0
set LLPXPGnum 0

set ULPOPEnum 43
set ULLPEAnum 1
set ULLPEBnum 2
set ULPOPGnum 18
set ULCLnum 0
set ULPXPEnum 0
set ULPXPGnum 0

#set X P
foreach X {A B D E H L P Y} {

	set mPE [mol new packmol/POPE.psf]
	mol addfile packmol/POPE.pdb
	set mLPEA [mol new packmol/LPE18.psf]
	mol addfile packmol/LPE18.pdb
	set mLPEB [mol new packmol/LPE16.psf]
	mol addfile packmol/LPE16.pdb
	set mPG [mol new packmol/POPG.psf]
	mol addfile packmol/POPG.pdb
	set mCL [mol new packmol/CL.psf]
	mol addfile packmol/CL.pdb


	if {"$X" != "P"} {

		set mPXE [mol new packmol/PEPGs/P${X}PE.psf]
		mol addfile packmol/PEPGs/P${X}PE.pdb
		set mPXG [mol new packmol/PEPGs/P${X}PG.psf]
		mol addfile packmol/PEPGs/P${X}PG.pdb

		set LLPOPEnum 27
		set LLPXPEnum 28
		set LLPOPGnum 2
		set LLPXPGnum 2
		set ULPOPEnum 21
		set ULPXPEnum 22
		set ULPOPGnum 9
		set ULPXPGnum 9
	} else {
		set LLPOPEnum 55
		set LLLPEAnum 0
		set LLLPEBnum 1
		set LLPOPGnum 4
		set LLCLnum 2
		set LLPXPEnum 0
		set LLPXPGnum 0
	
		set ULPOPEnum 43
		set ULLPEAnum 1
		set ULLPEBnum 2
		set ULPOPGnum 18
		set ULCLnum 0
		set ULPXPEnum 0
		set ULPXPGnum 0
	}


	set IM [mol new packmol/results/IM64${X}2.pdb]
	set sel [atomselect $IM "resname LLE"]
	$sel set resname CL
	set sel [atomselect $IM "resname LLA"]
	

	if {1} {
	###set PDB segnames for psf files...
	##
	## (sometimes) for LPEA or CL residue must be resid because packmol drops bonds
	## 	
	#### LOWER LEAFLET
		set sel [atomselect $IM "resname LPEA and same resid as (name P and z < 0)"]
		$sel set segname LLA
		set sel [atomselect $IM "resname LPEB and same residue as (name P and z < 0)"]
		$sel set segname LLB
		set sel [atomselect $IM "resname POPE and same residue as (name P and z < 0)"]
		$sel set segname LLC
		set sel [atomselect $IM "resname POPG and same residue as (name P and z < 0)"]
		$sel set segname LLD
		set sel [atomselect $IM "resname CL and same resid as (name P1 and z < 0)"]
		$sel set segname LLE
	
		set sel [atomselect $IM "resname PAPE PBPE PDPE PEPE PHPE PLPE PYPE and same resid as (name P and z < 0)"]
	$sel set segname LLF
	
	set sel [atomselect $IM "resname PAPG PBPG PDPG PEPG PHPG PLPG PYPG and same resid as (name P and z < 0)"]
	$sel set segname LLG
	
	#### UPPER LEAFLET
		set sel [atomselect $IM "resname LPEA and same resid as (name P and z > 0)"]
		$sel set segname ULA
		set sel [atomselect $IM "resname LPEB and same residue as (name P and z > 0)"]
		$sel set segname ULB
		set sel [atomselect $IM "resname POPE and same residue as (name P and z > 0)"]
		$sel set segname ULC
		$sel writepdb test.pdb
		set sel [atomselect $IM "resname POPG and same residue as (name P and z > 0)"]
		$sel set segname ULD
		set sel [atomselect $IM "resname CL and same residue as (name P1 and z > 0)"]
		$sel set segname ULE
	
		set sel [atomselect $IM "resname PAPE PBPE PDPE PEPE PHPE PLPE PYPE and same residue as (name P and z > 0)"]
		$sel set segname ULF
		
		set sel [atomselect $IM "resname PAPG PBPG PDPG PEPG PHPG PLPG PYPG and same residue as (name P and z > 0)"]
		$sel set segname ULG
		
		set all [atomselect $IM all] 
		$all writepdb IM64.pdb
	
	}

	if {1} {
		#####LOWER LEAFLET

		resetpsf
		set sel [atomselect $mLPEA all]
		$sel set segname LLA
		for {set i 1} {$i <= $LLLPEAnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf LPEA_only_LL.psf
	### LPEB (C:16) ####
		resetpsf
		set sel [atomselect $mLPEB all]
		$sel set segname LLB
		for {set i 1} {$i <= $LLLPEBnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf LPEB_only_LL.psf
	
	### POPE ####
		resetpsf
		set sel [atomselect $mPE all]
		$sel set segname LLC
		for {set i 1} {$i <= $LLPOPEnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf POPE_only_LL.psf
	
	### POPG ####
		resetpsf
		set sel [atomselect $mPG all]
		$sel set segname LLD
		for {set i 1} {$i <= $LLPOPGnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf POPG_only_LL.psf
	
	### CL ####
		resetpsf
		set sel [atomselect $mCL all]
		$sel set segname LLE
		for {set i 1} {$i <= $LLCLnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf CL_only_LL.psf
		
		if {"$X" != "P"} {	
	### PXPE ####
			resetpsf
			set sel [atomselect $mPXE all]
			$sel set segname LLF
			for {set i 1} {$i <= $LLPXPEnum} {incr i} {
				$sel set resid $i
				$sel writepsf tmp.psf
				$sel writepdb tmp.pdb
				readpsf tmp.psf
				coordpdb tmp.pdb
			}
			writepsf PXPE_only_LL.psf
			writepdb PXPE_only_LL.pdb
	
	
	### PXPG ####
			resetpsf
			set sel [atomselect $mPXG all]
			$sel set segname LLG
			for {set i 1} {$i <= $LLPXPGnum} {incr i} {
				$sel set resid $i
				$sel writepsf tmp.psf
				$sel writepdb tmp.pdb
				readpsf tmp.psf
				coordpdb tmp.pdb
			}
			writepsf PXPG_only_LL.psf
		}
	##############################################################
	##############################################################
	# UPPER LEAFLET (UL)
		
	### generate psf for each section
	
	### LPEA (C:18) ####
		resetpsf
		set sel [atomselect $mLPEA all]
		$sel set segname ULA
		for {set i 1} {$i <= $ULLPEAnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf LPEA_only_UL.psf
	### LPEB (C:16) ####
		resetpsf
		set sel [atomselect $mLPEB all]
		$sel set segname ULB
		for {set i 1} {$i <= $ULLPEBnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf LPEB_only_UL.psf
	
	### POPE ####
		resetpsf
		set sel [atomselect $mPE all]
		$sel set segname ULC
		for {set i 1} {$i <= $ULPOPEnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf POPE_only_UL.psf
	
	### POPG ####
		resetpsf
		set sel [atomselect $mPG all]
		$sel set segname ULD
		for {set i 1} {$i <= $ULPOPGnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf POPG_only_UL.psf
	
	### CL #### 
	#       >>> no CL in UL as is <<<
		resetpsf
		set sel [atomselect $mCL all]
		$sel set segname ULE
		for {set i 1} {$i <= $ULCLnum} {incr i} {
			$sel set resid $i
			$sel writepsf tmp.psf
			$sel writepdb tmp.pdb
			readpsf tmp.psf
			coordpdb tmp.pdb
		}
		writepsf CL_only_UL.psf
	
		if {"$X" != "P"} {
	### PXPE ####
			resetpsf
			set sel [atomselect $mPXE all]
			$sel set segname ULF
			for {set i 1} {$i <= $ULPXPEnum} {incr i} {
				$sel set resid $i
				$sel writepsf tmp.psf
				$sel writepdb tmp.pdb
				readpsf tmp.psf
				coordpdb tmp.pdb
			}
			writepsf PXPE_only_UL.psf
		
		
	### PXPG ####
			resetpsf
			set sel [atomselect $mPXG all]
			$sel set segname ULG
			for {set i 1} {$i <= $ULPXPGnum} {incr i} {
				$sel set resid $i
				$sel writepsf tmp.psf
				$sel writepdb tmp.pdb
				readpsf tmp.psf
				coordpdb tmp.pdb
			}
			writepsf PXPG_only_UL.psf
		}	
	
	}
	
	if {0} {
	######## PME
		set mPME [mol new peptides/PME/PME.psf]
		mol addfile peptides/PME/PME.pdb
		set selPME [atomselect $mPME all]
		$selPME moveby { 0 0 38 }
		
		rot all y 70
		rot all x 30
	
		$selPME writepdb PME.pdb
		$selPME writepsf PME.psf
	}
	
	if {1} {
		resetpsf
		foreach i {LPEA LPEB POPE POPG CL PXPE PXPG} {
			foreach j {LL UL} {
				if { [set ${j}${i}num] > 0} {
					readpsf ${i}_only_${j}.psf
				}
			}
		}
		#readpsf PME.psf
		coordpdb IM64.pdb
		#coordpdb PME.pdb
	
	
		writepdb IM64r.pdb
		writepsf IM64r.psf
	
		mol delete all
		mol new IM64r.psf
		mol addfile IM64r.pdb
	
	}
	
	if {1} {
	#### SOLVATING
		set IMsel [atomselect top "chain L U"]
	
		set r [measure minmax $IMsel]
		foreach {min max} $r { break }
		foreach {xmin ymin zmin} $min { break }
		foreach {xmax ymax zmax} $max { break }
		#set wbtop [subst {{ $xmin $ymin [expr $zmax - 5 ]} { $xmax $ymax [expr $zmax + 35]}}
	
	 
		set wbtop [subst {{ -35 -35 [expr $zmax - 2 ]} { 35 35 [expr $zmax + 35]}}]
		set wbbot [subst {{ -35 -35 [expr $zmin - 35 ]} { 35 35 [expr $zmin + 2]}}]
		solvate IM64r.psf IM64r.pdb -minmax $wbtop -b 2 -o tmpt
		solvate tmpt.psf tmpt.pdb -minmax $wbbot -s WA -b 2 -o IM64_wb
	
	}
	if {1} {
		autoionize -psf IM64_wb.psf -pdb  IM64_wb.pdb -sc 0.15 -cation POT -o IM${X}_wbi
	}		

	if {1} {
###########################################
####		no SMD
####
		
		exec sed -i "s/XXX/${X}/g" PMEPUFA/${X}b/s0.inp
		exec sed -i "s/XXX/${X}/g" PMEPUFA/${X}b/s0.qsub

 		for {set i 1} {$i <= 10} {incr i} {
			exec sed "s/restart 0/restart ${i}/g" PMEPUFA/${X}b/s0.inp > PMEPUFA/${X}b/s${i}.inp
			exec sed "s/s0/s${i}/g" PMEPUFA/${X}b/s0.qsub > PMEPUFA/${X}b/s${i}.qsub
			exec sed -i "s/b.0/b.${i}/g" PMEPUFA/${X}b/s${i}.qsub
		}			
		
		
		

		mol new IM${X}_wbi.psf
		mol addfile IM${X}_wbi.pdb waitfor all
	
		set sel [atomselect top all]
		set beta 0
		set occupancy 0

		$sel writepsf PMEPUFA/${X}b/${X}.psf
		$sel writepdb PMEPUFA/${X}b/${X}.pdb
		 		
	}
}


mol delete all

foreach X {A B D E H L P Y} {
	mol new PMEPUFA/${X}b/${X}.psf
	mol addfile PMEPUFA/${X}b/${X}.pdb
}
