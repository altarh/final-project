echo "making 324976703_206667248_212875066_project dir..."
mkdir 324976703_206667248_212875066_project
echo "copying files into dir:"
echo "symnmf.py..."
cp symnmf.py 324976703_206667248_212875066_project
echo "symnmf.c..."
cp symnmf.c 324976703_206667248_212875066_project
echo "symnmfmodule.c..."
cp symnmfmodule.c 324976703_206667248_212875066_project
echo "symnmf.h..."
cp symnmf.h 324976703_206667248_212875066_project
echo "analysis.py..."
cp analysis.py 324976703_206667248_212875066_project
echo "setup.py..."
cp setup.py 324976703_206667248_212875066_project
echo "Makefile..."
cp Makefile 324976703_206667248_212875066_project
echo "zipping..."
tar -czvf 324976703_206667248_212875066_project.tar.gz 324976703_206667248_212875066_project
echo "deleting dir..."
rm -r 324976703_206667248_212875066_project