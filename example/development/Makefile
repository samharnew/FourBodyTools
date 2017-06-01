
include $(PWD)/Makefile.arch

HdrSuf = h
HdrSufOther = hpp
SrcSuf = cpp
ToyDir = $(PWD)/../..
MintDir      = $(HOME)/Documents/D0MixingAnalysis/sublime_local/MINT2_dargent/Mint2/src/Mint
MojitoDir    = $(HOME)/Documents/D0MixingAnalysis/sublime_local/MINT2_dargent/Mint2/src/Mojito
HyperPlotDir = $(HOME)/Documents/Work/Code/HyperPlot


MOREINCS := $(ToyDir)/include
MOREINCS += $(wildcard $(MintDir)/*/src)
MOREINCS += $(wildcard $(MojitoDir)/*/src)
MOREINCS += $(wildcard $(MojitoDir)/../../)
MOREINCS += $(wildcard $(HyperPlotDir)/include)

MOREINCSFLAGS := $(patsubst %,-I%,$(MOREINCS))

CXXFLAGS += $(MOREINCSFLAGS)
DEPCXXFLAGS := CXXFLAGS


LIBS += -lMinuit -lRooFit -lRooStats -lRooFitCore -lHtml 


#pipe file names containing main script to sed to remove extenstion
PROG = $(shell grep -l main src/*.cpp | sed 's/\.cpp//g')

SRCS := $(wildcard $(ToyDir)/src/*.$(SrcSuf))      #create list of .cc files


SRCS += $(wildcard $(PWD)/src/*.$(SrcSuf))
HDRS := $(wildcard $(ToyDir)/include/*.$(HdrSuf))  #create list of .h  files

SRCS += $(wildcard $(MintDir)/*/src/*.$(SrcSuf))
SRCS += $(wildcard $(MojitoDir)/*/src/*.$(SrcSuf))
HDRS += $(wildcard $(MintDir)/*/src/*.$(HdrSuf))
HDRS += $(wildcard $(MintDir)/../../Mint/*.$(HdrSuf))
HDRS += $(wildcard $(MojitoDir)/*/src/*.$(HdrSuf))

SRCS += $(wildcard $(HyperPlotDir)/src/*.$(SrcSuf))      #create list of .cc files
HDRS _= $(wildcard $(HyperPlotDir)/include/*.$(HdrSuf))  #create list of .h  files


OBJS := $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(SRCS))  #create a list of .o files from .cc files


SRCSLOCAL += $(wildcard src/*.$(SrcSuf))

OBJSLOCAL := $(patsubst %.$(SrcSuf),%.$(ObjSuf),$(SRCSLOCAL))


# top-level rule, to compile everything.
all: $(PROG)
	@echo ""; echo " ============== make all done ==============="; echo ""

# don't link, just compile
src: $(OBJSLOCAL)
	@echo ""; echo " ============== make src done ==============="; echo ""

depend: .depend
	@echo ""; echo " ============== make depend done ============"; echo ""

.depend: $(SRCS) $(HDRS)
	@echo ""; echo " ===== Figuring out dependencies ============"; echo ""
	@makedepend -f- -- $(DEPCXXFLAGS) -- $(SRCS) > .depend 2> .makedepend.err
	@echo ""; echo " result in .depend, errors in .makedepend.err"; echo ""
	@echo ""; echo " ======== Dependencies done. ================"; echo ""

clean: 
	$(RM) $(OBJS) $(PROG) .depend .makedepend.err
	@echo ""; echo " =============== make clean done ============"; echo ""

debug:
	@echo "sources " $(SRCS)
	@echo "objects " $(OBJS)
	@echo "pwd " $(PWD)
	@echo "CXXFLAGS " $(CXXFLAGS)
	@echo "MintDir " $(MintDir)
	
$(PROG): $(OBJS)
	@echo ""; echo " =============== linking ===================="; echo ""
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $(PROG)
	@echo ""; echo " =============== linking done ==============="; echo ""

include .depend

%.(ObjSuf): %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<
##
# DO NOT DELETE
