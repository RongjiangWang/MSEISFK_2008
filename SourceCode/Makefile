FORFILES = \
caxcb.f      mslayer.f    msmatrix.f   timeprt.f     cgemp.f      mswaveno.f  \
mssublay.f   mskern.f     mshkwa.f     msfkmain.f    cmemcpy.f    gemp.f      \
mshkpsv.f    msmatinv.f   msqmodel.f   getdata.f

OBJECTS =	$(FORFILES:%.f=%.o)

FFLAGS =	-fast -O3 -silent

PROGRAM =	mseisfk

$(PROGRAM): 	$(OBJECTS)
		$(FC) $(FFLAGS) $(OBJECTS) -o $@
caxcb.o:caxcb.f
mslayer.o:msglobal.h mslayer.f
msmatrix.o:msglobal.h msmatrix.f
timeprt.o:timeprt.f
cgemp.o:cgemp.f
mswaveno.o:msglobal.h mswaveno.f
msmatinv.o:msglobal.h msmatinv.f
mssublay.o:msglobal.h mssublay.f
msfkmain.o:msglobal.h msfkmain.f
msqmodel.o:msglobal.h msqmodel.f
cmemcpy.o:cmemcpy.f
gemp.o:gemp.f
mskern.o:msglobal.h mskern.f
mshkpsv.o:msglobal.h mshkpsv.f
mshkwa.o:msglobal.h mshkwa.f
getdata.o:getdata.f









