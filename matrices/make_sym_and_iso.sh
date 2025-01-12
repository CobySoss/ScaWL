# Check if two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory> <filename>"
    exit 1
fi

./gen_symmetric ./$1/$2.mtx ./$1/$2_sym.mtx
./gen_iso_relabel ./$1/$2_sym.mtx ./$1/$2_sym_iso.mtx