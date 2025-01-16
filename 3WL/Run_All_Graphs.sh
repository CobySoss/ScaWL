#!/bin/sh
sh run_tests_entire_set.sh ../matrices/494_bus 494_bus_sym.mtx  ../matrices/494_bus 494_bus_sym_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/662_bus 662_bus_sym.mtx  ../matrices/662_bus 662_bus_sym_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/685_bus 685_bus_sym.mtx  ../matrices/685_bus 685_bus_sym_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_multinode_driver.sh ../matrices/b_dyn b_dyn.mtx ../matrices/b_dyn b_dyn_iso.mtx '1-12:00:00'
sh run_tests_multinode_driver.sh ../matrices/can_1072 can_1072_sym.mtx ../matrices/can_1072 can_1072_sym_iso.mtx '1-12:00:00'
sh run_tests_multinode_driver.sh ../matrices/msc01050 msc01050_sym.mtx ../matrices/msc01050 msc01050_sym_iso.mtx '1-12:00:00'
sh run_tests_multinode_driver.sh ../matrices/1138_bus 1138_bus_sym.mtx ../matrices/1138_bus 1138_bus_sym_iso.mtx '1-12:00:00'
sh run_tests_multinode_driver.sh ../matrices/bcsstk08 bcsstk08_sym.mtx ../matrices/bcsstk08 bcsstk08_sym_iso.mtx '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/bcsstk27 bcsstk27_sym.mtx ../matrices/bcsstk27 bcsstk27_sym_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/celegansneural celegansneural.mtx ../matrices/celegansneural celegansneural_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/fs_680_3 fs_680_3.mtx ../matrices/fs_680_3 fs_680_3_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/fs_760_1 fs_760_1.mtx ../matrices/fs_760_1 fs_760_1_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/gr_30_30 gr_30_30_sym.mtx ../matrices/gr_30_30 gr_30_30_sym_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/rbsa480 rbsa480.mtx ../matrices/rbsa480 rbsa480_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/G12 G12_sym.mtx ../matrices/G12 G12_sym_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/ex22 ex22.mtx ../matrices/ex22 ex22_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/qh882 qh882.mtx ../matrices/qh882 qh882_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/orsirr_2 orsirr_2.mtx ../matrices/orsirr_2 orsirr_2_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/bp_1000 bp_1000.mtx ../matrices/bp_1000 bp_1000_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/saylr3 saylr3_sym.mtx ../matrices/saylr3 saylr3_sym_iso.mtx '1-12:00:00' '1-12:00:00'
sh run_tests_entire_set.sh ../matrices/cdde3 cdde3.mtx ../matrices/cdde3 cdde3_iso.mtx '1-12:00:00' '1-12:00:00'