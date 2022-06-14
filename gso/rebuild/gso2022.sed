s/REAL    /DOUBLE PRECISION/g
s/PROC,//g
s/^[ \t]*EXTERNAL/C     EXTERNAL/
/CALL PROC/s/^./C/g
s/^\t/      /g
