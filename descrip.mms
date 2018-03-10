cflags = $(cflags)/optimize=inline
rtlib = sys$library:vaxcrtl

all : procoseq.exe rcs.exe
	@ write sys$output "Done."

procoseq.exe : procoseq.obj fperiod.obj longer.obj inumber.obj
        $(link)$(linkflags)/notrace procoseq.obj, fperiod.obj, longer.obj, inumber.obj, $(rtlib)/lib

rcs.exe : rcs.obj fperiod.obj longer.obj inumber.obj
        $(link)$(linkflags)/notrace rcs.obj, fperiod.obj, longer.obj, inumber.obj, $(rtlib)/lib

abc.exe : abc.obj longer.obj inumber.obj
        $(link)$(linkflags)/notrace abc.obj, longer.obj, inumber.obj, $(rtlib)/lib

.c.obj :
        $(cc) $(cflags) $(mms$source)
