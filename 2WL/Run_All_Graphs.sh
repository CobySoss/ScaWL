#!/bin/sh
sh run_tests_entire_set.sh ../matrices/appu appu.mtx ../matrices/appu appu_iso.mtx 
sh run_tests_entire_set.sh ../matrices/bcsstk25 bcsstk25_sym.mtx ../matrices/bcsstk25 bcsstk25_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/bcsstk36 bcsstk36_sym.mtx ../matrices/bcsstk36 bcsstk36_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/bcsstk32 bcsstk32_sym.mtx ../matrices/bcsstk32 bcsstk32_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/coater2 coater2.mtx ../matrices/coater2 coater2_iso.mtx
sh run_tests_entire_set.sh ../matrices/dw8192 dw8192.mtx ../matrices/dw8192 dw8192_iso.mtx
sh run_tests_entire_set.sh ../matrices/epb1 epb1.mtx ../matrices/epb1 epb1_iso.mtx
sh run_tests_entire_set.sh ../matrices/G61 G61_sym.mtx ../matrices/G61 G61_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/hangGlider_4 hangGlider_4_sym.mtx ../matrices/hangGlider_4 hangGlider_4_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/Kaufhold Kaufhold.mtx ../matrices/Kaufhold Kaufhold_iso.mtx
sh run_tests_multinode_driver.sh ../matrices/nd12k nd12k_sym.mtx ../matrices/nd12k nd12k_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/olafu olafu_sym.mtx ../matrices/olafu olafu_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/poli_large poli_large.mtx ../matrices/poli_large poli_large_iso.mtx
sh run_tests_multinode_driver.sh ../matrices/ct20stif ct20stif_sym.mtx ../matrices/ct20stif ct20stif_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/bayer04 bayer04.mtx ../matrices/bayer04 bayer04_iso.mtx
sh run_tests_multinode_driver.sh ../matrices/bcsstm35 bcsstm35_sym.mtx ../matrices/bcsstm35 bcsstm35_sym_iso.mtx
sh run_tests_multinode_driver.sh ../matrices/onetone1 onetone1.mtx ../matrices/onetone1 onetone1_iso.mtx
sh run_tests_entire_set.sh ../matrices/crystk03 crystk03_sym.mtx ../matrices/crystk03 crystk03_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/g7jac080sc g7jac080sc.mtx ../matrices/g7jac080sc g7jac080sc_iso.mtx
sh run_tests_multinode_driver.sh ../matrices/lhr34 lhr34.mtx ../matrices/lhr34 lhr34_iso.mtx
sh run_tests_multinode_driver.sh ../matrices/pwt pwt_sym.mtx ../matrices/pwt pwt_sym_iso.mtx
sh run_tests_entire_set.sh ../matrices/raefsky3 raefsky3.mtx ../matrices/raefsky3 raefsky3_iso.mtx
sh run_tests_multinode_driver.sh ../matrices/wang3 wang3.mtx ../matrices/wang3 wang3_iso.mtx
sh run_tests_entire_set.sh ../matrices/pkustk01 pkustk01_sym.mtx ../matrices/pkustk01 pkustk01_sym_iso.mtx

