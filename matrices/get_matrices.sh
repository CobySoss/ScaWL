if [ -n "$1" ] && [ "$1" != "no_cert_validation" ]; then
  echo "Invalid parameter"
  exit 1
fi

if [ "$1" == "no_cert_validation" ]; then
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Brethour/coater2.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Simon/appu.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk32.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Boeing/bcsstk36.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Boeing/ct20stif.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Bai/dw8192.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Averous/epb1.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk25.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/VDOL/hangGlider_4.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/MathWorks/Kaufhold.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/ND/nd12k.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Simon/olafu.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Grund/poli_large.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Gset/G61.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Grund/b_dyn.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/1138_bus.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk08.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/494_bus.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/662_bus.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/685_bus.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Boeing/msc01050.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/can_1072.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk27.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Newman/celegansneural.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/fs_680_3.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/fs_760_1.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/gr_30_30.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Bai/rbsa480.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Gset/G12.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/FIDAP/ex22.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Bai/qh882.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/orsirr_2.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/bp_1000.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/HB/saylr3.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Bai/cdde3.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Grund/bayer04.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Boeing/bcsstm35.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/ATandT/onetone1.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Boeing/crystk03.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Hollinger/g7jac080sc.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Mallya/lhr34.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Nasa/pwt.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Simon/raefsky3.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Wang/wang3.tar.gz
  wget --no-check-certificate https://www.cise.ufl.edu/research/sparse/MM/Chen/pkustk01.tar.gz    
else
  wget https://www.cise.ufl.edu/research/sparse/MM/Brethour/coater2.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Simon/appu.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk32.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Boeing/bcsstk36.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Boeing/ct20stif.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Bai/dw8192.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Averous/epb1.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk25.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/VDOL/hangGlider_4.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/MathWorks/Kaufhold.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/ND/nd12k.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Simon/olafu.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Grund/poli_large.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Gset/G61.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Grund/b_dyn.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/1138_bus.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk08.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/494_bus.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/662_bus.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/685_bus.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Boeing/msc01050.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/can_1072.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/bcsstk27.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Newman/celegansneural.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/fs_680_3.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/fs_760_1.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/gr_30_30.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Bai/rbsa480.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Gset/G12.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/FIDAP/ex22.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Bai/qh882.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/orsirr_2.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/bp_1000.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/HB/saylr3.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Bai/cdde3.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Grund/bayer04.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Boeing/bcsstm35.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/ATandT/onetone1.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Boeing/crystk03.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Hollinger/g7jac080sc.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Mallya/lhr34.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Nasa/pwt.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Simon/raefsky3.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Wang/wang3.tar.gz
  wget https://www.cise.ufl.edu/research/sparse/MM/Chen/pkustk01.tar.gz    
fi

for file in *.gz; do
  if [ -f "$file" ]; then
    echo "Untarring $file..."
    tar -xvf "$file"
  fi
done
g++ -std=c++11 ../util/bidirection_adder.c -o gen_symmetric
g++ -std=c++11 ../util/iso_create.c -o gen_iso_relabel

./make_sym_and_iso.sh 494_bus 494_bus
./make_sym_and_iso.sh 662_bus 662_bus
./make_sym_and_iso.sh 685_bus 685_bus
./make_sym_and_iso.sh 1138_bus 1138_bus
./make_iso.sh b_dyn b_dyn
./make_sym_and_iso.sh can_1072 can_1072
./make_iso.sh celegansneural celegansneural
./make_iso.sh coater2 coater2
./make_sym_and_iso.sh ct20stif ct20stif
./make_iso.sh dw8192 dw8192
./make_iso.sh epb1 epb1
./make_iso.sh fs_680_3 fs_680_3
./make_iso.sh fs_760_1 fs_760_1
./make_sym_and_iso.sh gr_30_30 gr_30_30
./make_sym_and_iso.sh hangGlider_4 hangGlider_4
./make_iso.sh Kaufhold Kaufhold
./make_sym_and_iso.sh msc01050 msc01050
./make_sym_and_iso.sh nd12k nd12k
./make_sym_and_iso.sh olafu olafu
./make_iso.sh poli_large poli_large
./make_iso.sh rbsa480 rbsa480
./make_iso.sh appu appu
./make_sym_and_iso.sh bcsstk08 bcsstk08
./make_sym_and_iso.sh bcsstk25 bcsstk25
./make_sym_and_iso.sh bcsstk27 bcsstk27
./make_sym_and_iso.sh bcsstk32 bcsstk32
./make_sym_and_iso.sh bcsstk36 bcsstk36
./make_sym_and_iso.sh G61 G61
./make_sym_and_iso.sh G12 G12
./make_iso.sh ex22 ex22
./make_iso.sh qh882 qh882
./make_iso.sh orsirr_2 orsirr_2
./make_iso.sh bp_1000 bp_1000
./make_sym_and_iso.sh saylr3 saylr3
./make_iso.sh cdde3 cdde3 
./make_iso.sh bayer04 bayer04
./make_sym_and_iso.sh bcsstm35 bcsstm35
./make_iso.sh onetone1 onetone1
./make_sym_and_iso.sh crystk03 crystk03
./make_iso.sh g7jac080sc g7jac080sc
./make_iso.sh lhr34 lhr34
./make_sym_and_iso.sh pwt pwt
./make_iso.sh raefsky3 raefsky3
./make_iso.sh wang3 wang3
./make_sym_and_iso.sh pkustk01 pkustk01

echo "Finished obtaining and generating matrices"

