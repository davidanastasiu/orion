rm -rf $1
svn up
svn export http://dminers.dtc.umn.edu/svn/programs/david/trunk/orion/trunk $1
svn export http://dminers.dtc.umn.edu/svn/programs/david/trunk/orion/data $1/data
mkdir $1/manual
svn export http://dminers.dtc.umn.edu/svn/programs/david/trunk/orion/manual/trunk/out/manual.pdf $1/manual/manual.pdf
svn export http://dminers.dtc.umn.edu/svn/programs/david/trunk/orion/manual/Orion_ICDE15.pdf $1/manual/Orion_ICDE15.pdf
svn export http://dminers.dtc.umn.edu/svn/libs/GKlib/trunk $1/GKlib
rm -rf $1/TODO
rm -rf $1/test
rm -rf $1/utils
tar -czf $1.tar.gz $1
rm -rf $1
