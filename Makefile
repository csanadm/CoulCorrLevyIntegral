PROGRAMS = $(basename $(shell ls *.cc))

WrkDir := $(shell pwd)
SrcDir        = .
ObjDir        = .
DepDir        = ./deps
ExeDir        = .

CXX           = g++
CFLAGS        = -O3 -Wall -fPIC -fno-inline -std=c++14
LD            = g++
LDFLAGS       = -O3 

COMMON_SOURCES = $(shell ls *.cpp)
COMMON_OBJECTS = $(addprefix $(ObjDir)/, $(addsuffix .o,$(notdir $(basename $(COMMON_SOURCES)))))
SOURCES = $(addprefix $(SrcDir)/,$(addsuffix .cc,$(PROGRAMS))) $(COMMON_SOURCES)
ALL_SOURCES = $(sort $(SOURCES))

define COMPILE_TEMPLATE
-include $(DepDir)/$(notdir $(basename $(1))).d
$(ObjDir)/$(notdir $(basename $(1))).o:
	@echo ""
	$(CXX) $(CFLAGS) -c -MD -MP -MF $(DepDir)/$(notdir $(basename $(1))).d $(1) -o $$@
endef
$(foreach source, $(ALL_SOURCES), $(eval $(call COMPILE_TEMPLATE,$(source))))

define FINAL_LINK_TEMPLATE
$(ExeDir)/$(1).exe: $(2) $(ObjDir)/$(1).o
	@echo ""
	$(LD) $(LDFLAGS) $$^ $(SYSLIBS) -o $$@
endef
$(foreach program, $(PROGRAMS), $(eval $(call FINAL_LINK_TEMPLATE,$(program),$(COMMON_OBJECTS))))

all: $(addprefix $(ExeDir)/,$(addsuffix .exe, $(PROGRAMS)))

clean:
	@rm -f $(ExeDir)/*.exe $(ObjDir)/*.o $(DepDir)/*.d
