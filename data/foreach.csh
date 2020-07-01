#/bin/tcsh
cd /x/zhangc/tmp/a/

rm a.lis
foreach fl (*.ent)
edpdb $fl << eof > a.jnk
file
residu arg | side | atom except h* 1* 2*; gr a
dist a 2.0 3.3 1
quit
eof
grep -i "\.ent" a.jnk >> a.lis
grep -i ",d=" a.jnk >> a.lis
end
