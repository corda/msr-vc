# Switch here between the two matching curves used by Geppetto in
# normal (inner) and bootstrapping (outer) mode.
#
# The "tiny" curves use a 54-bit prime instead of a 254-bit prime; 
# their processing is faster but insecure; they are just for testing.
#
TINY=
# TINY=TINY

ifeq ($(TINY),TINY)
INNER=7  # Tiny BN
OUTER=8  # Tiny CP
else
INNER=3  # Secure BN
OUTER=5  # Secure CP
endif

# Switch here between the Release and Debug builds of Geppetto;
# be sure to build the corresponding configuration of the VS solution.
# (Debug helps debugging Geppetto C++ libraries but is much slower.)
#
BLT=Release
# BLT=Debug

# ------------- external dependencies

# Using the clang C compiler (tested with LLVM 3.5.0, 3.5.2, 3.6.2)
# clang/LLVM should be installed independently, 
# as the cygwin packages generate old unsupported LLVM syntax
clang = /cygdrive/c/Program\ Files\ \(x86\)/LLVM/bin/clang -O2

# for compiling llvm input to Geppetto
clangs = $(clang) -S -emit-llvm -DMQAP=1 -fno-slp-vectorize -fno-vectorize 

# for compiling to native code (and comparing results)
clangw = $(clang) -D_WIN32 -DNATIVE -I"C:/Cygwin/usr/include" '-D__warning__(arg)= '

# Using fslex and fsyacc for LLVM syntax
fs4bin = /cygdrive/c/Program\ Files\ \(x86\)/FSharpPowerPack-4.0.0.0/bin
build = ./build

# ------------- Geppetto 

PERF=
core = ../bin/$(BLT)/geppetto.exe --fflib v 0 $(PERF) $(OPT)
base = $(core) --encoding $(INNER)
outer = $(core) --encoding $(OUTER)

# ------------- QAP-friendly cryptographic libraries (for bootstrapping)

# Qapped Arith, for native and geppetto compilation
ARITH=../src/arith/arith/
QAP_ARITH=../src/arith-qap/arith/Qapped/

# ------------- main targets for testing

all:: test-all cli-sample wp-sample bootstrap-sample

cli-sample:: fg matrixes-small matrixes-one bankFile branch sha-all

wp-sample:: fgw matrix.wp stages

bootstrap-sample:: whole.wp whole.verify wide-loop.wp wide-loop.verify

test-all:: tests trivials 

matrixes-one: one-matrix shared overlap


build: 
	mkdir -p $(build)

build-clean:
	rm -rf build/
	mkdir build

clean: build-clean

# ------------- building from VS

msbuild = "/cygdrive/c/Program Files (x86)/MSBuild/12.0/Bin/MSBuild.exe"
# vcbuild = $(msbuild) /p:Configuration=$(1) ../windows/compiler.sln
# vcclean = $(msbuild) /v:minimal /p:Configuration=$(1) /t:Clean ../windows/compiler.sln

SRC=../src/fsharp

# Generating F# source files

$(SRC)/QAPEval.fs: $(SRC)/QAPGen.fs
	cd $(SRC); cpp QAPGen.fs | grep -v '^#' > QAPEval.fs

llvm:: $(SRC)/LLVM_grammar.fs $(SRC)/LLVM_syntax.fs $(SRC)/QAPEval.fs
	@echo "*** recompile in VS ***"

# now automated in VS:
# $(SRC)/LLVM_grammar.fs $(SRC)/LLVM_syntax.fs: $(SRC)/LLVM_grammar.fsy $(SRC)/LLVM_syntax.fsl
#	cd $(SRC) ; $(fs4bin)/fsyacc.exe LLVM_grammar.fsy --module LLVM_grammar
#	cd $(SRC) ; $(fs4bin)/fslex.exe LLVM_syntax.fsl --unicode

msbuild: build ../windows/compiler.sln llvm
	$(msbuild) /p:Configuration=Release ../windows/compiler.sln
	$(msbuild) /p:Configuration=Debug ../windows/compiler.sln

msclean: clean
	$(msbuild) /p:Configuration=Release /t:Clean ../windows/compiler.sln
	$(msbuild) /p:Configuration=Debug /t:Clean ../windows/compiler.sln

# ------------- basic tests of FFLib/F# interop 

tests:: field.test fg.test simple.test split.test

%.test:: 
	$(core) $(ENCODING) -v --test $*

fields::
	(for e in 1 3 4 5 6 7 8; do $(core) -v --encoding $$e ; done) | grep setting

#--------------------------------------------------------------------
# default entries for running Geppetto in release or debug mode, in
# command-line style, assuming fixed outsource(Input,Output)

%.run : %.s
	$(base) --load $*.s --compile outsource --keygen --input Input random --eval outsource 

%.dbg : %.s
	$(base) -v --check --load $*.s --compile outsource --keygen --input Input random --eval outsource --output Input --output Output

$(build)/%.0.c : %.c build
	cpp -DMQAP=1 -DPARAM=0 $< -o $@
$(build)/%.1.c : %.c build
	cpp -DMQAP=1 -DPARAM=1 $< -o $@
$(build)/%.2.c : %.c build
	cpp -DMQAP=1 -DPARAM=2 $< -o $@
$(build)/%.3.c : %.c build
	cpp -DMQAP=1 -DPARAM=3 $< -o $@
$(build)/%.4.c : %.c build
	cpp -DMQAP=1 -DPARAM=4 $< -o $@
$(build)/%.5.c : %.c build
	cpp -DMQAP=1 -DPARAM=5 $< -o $@
$(build)/%.6.c : %.c build
	cpp -DMQAP=1 -DPARAM=6 $< -o $@
$(build)/%.7.c : %.c build
	cpp -DMQAP=1 -DPARAM=7 $< -o $@
$(build)/%.8.c : %.c build
	cpp -DMQAP=1 -DPARAM=8 $< -o $@
$(build)/%.99.c : %.c build
	cpp -DMQAP=1 -DPARAM=99 $< -o $@

$(build)/%.s : $(build)/%.c
	$(clangs) $<  -o $@

$(build)/%.s : %.c
	$(clangs) $<  -o $@

# ---------------------------------------------
# simple one-matrix and two-matrix examples

build/one-matrix.exe: one-matrix-ifc.h one-matrix.c
	$(clangw) -DPARAM=0 -DTRUNCATE=1 one-matrix.c -o $@

one-matrix: build/one-matrix.0.s build/one-matrix.exe
	build/one-matrix.exe
	$(base) -q \
--load build/one-matrix.0.s --compile outsource 1,2,3,4 --keygen \
--input A 5,6,7,8 --eval outsource 1,2,3,4 --output R 

matrixes-small:: \
build/two-matrix.0.run \
build/two-matrix.1.run \
build/two-matrix.2.run \
build/two-matrix.3.run \
build/two-matrix.4.run 

matrixes-all:: \
build/two-matrix.0.run \
build/two-matrix.1.run \
build/two-matrix.2.run \
build/two-matrix.3.run \
build/two-matrix.4.run \
build/two-matrix.5.run 

matrixes-sparse:: \
build/two-matrix.2.run \
build/two-matrix.4.run \
build/two-matrix.5.run 

# ---------------------------------------------
# a basic 2-QAP example evaluating "g(f(2,3))"
# with explicit command lines to compile and run f & g.

fg:: build/mqap_fg.s
	$(base) -e --check $(ENCODING) \
--load build/mqap_fg.s \
--compile f 4 \
--compile g \
--keygen  \
--input input 2,3 --eval f 4 \
--eval g \
--output output

# whole-program variant of the same example. 
# the OPAQUE flag optionally makes the intermediate bank opaque
# (so that it is verified rather than recomputed by the verifier)

fgw:
	make OPT="-e" mqap_fgw.wp
	make OPT="-e" COPT="-DOPAQUE" mqap_fgw.wp

# ---------------------------------------------
# some other m-QAP examples
# 
shared::  build/one-matrix.0.s
	$(base) --check $(ENCODING) --load build/one-matrix.0.s \
--compile outsource 1,2,0,1 \
--compile outsource 2,1,1,0 \
--keygen --input A 3,5,7,11 --eval outsource 1,2,0,1 --output R

overlap::  build/one-matrix.0.s
	$(base)  $(ENCODING) --check --load build/one-matrix.0.s \
--compile outsource 1,2,0,1 \
--compile outsource 2,1,1,0 \
--keygen --input A random \
--eval outsource 1,2,0,1 --output R


# ---------------------------------------------------------------------------
# We implement SHA1 and SHA256, with variable-length inputs and
# padding.  This provides another easily-tunable benchmark, with a
# "dense" QAP whose values are almost all bits (hence compilation and
# polynomial division are relatively costly).
#
# The first tests are fixed length, single compression function.
# The others rely on a shared library (sha-library.h)
#
# For SHA1, we use mostly PARAM=1; (PARAM=0 is deprecated; PARAM=2 has
# different hand-tuned C code for the inner bitwise operations, which
# does not help.)

# for regression testing
sha-all:: sha.1 sha.2 sha1 sha256 shaboth
	make LENGTH=67 sha1random
	make LENGTH=67 sha256random

build/sha.1.exe build/sha.2.exe: pinocchio.h pinocchio.c sha-ifc.h sha.c 
	$(clangw) -D_WIN32 -DPARAM=1 sha.c  -o $@
	$(clangw) -D_WIN32 -DPARAM=2 sha.c  -o $@

# test vector from Pinocchio v1.
SHA_ABC = --input Input 1633837824,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 --output Input

sha.1:: build/sha.1.exe build/sha.1.s  
	build/sha.1.exe
	$(base) --load build/sha.1.s --compile outsource --keygen $(SHA_ABC) --eval outsource --output Input --output Output

# still sha1, with a hand-optimized C implementation
sha.2:: build/sha.2.exe build/sha.2.s
	build/sha.2.exe
	$(base) --load build/sha.2.s --compile outsource --keygen $(SHA_ABC) --eval outsource --output Input --output Output

# variable-length SHA1 and SHA256 with padding; 96 bytes require two blocks, 159 require 3.
LENGTH=3
#LENGTH=96
#LENGTH=159
SHA_INPUT_3 = 1633837824
SHA_INPUT_96 = 1633837824,0,0,0,0,0,0,0,0,345,0,0,0,23423,0,2423,0,0,0,0,0,0,5,42
SHA_INPUT_159 = 1633837824,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1633837824,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

sha1::  
	$(clangs) -DPARAM=1 -DLENGTH=$(LENGTH) sha.c -o build/sha.3.s 
	$(base) \
--load build/sha.3.s  --compile sha1n $(LENGTH) --keygen \
--input Msg $(SHA_INPUT_$(LENGTH)) --eval sha1n $(LENGTH) --output Output

sha256::  
	$(clangs) -DPARAM=1 -DLENGTH=$(LENGTH) sha.c -o build/sha.3.s
	$(base) \
--load build/sha.3.s  --compile sha256n $(LENGTH) --keygen \
--input Msg $(SHA_INPUT_$(LENGTH)) --eval sha256n $(LENGTH) --output Output256

# confirming the low overhead of MQAPs (notice a slight increase in the runtime for the smaller SHA1 QAP).
shaboth::
	$(clangs) -DPARAM=1 -DLENGTH=$(LENGTH) sha.c -o build/sha.3.s
	$(base) \
--load build/sha.3.s --compile sha1n $(LENGTH) --compile sha256n $(LENGTH) --keygen \
--input Msg $(SHA_INPUT_$(LENGTH)) --eval sha1n $(LENGTH) --eval sha256n $(LENGTH) 

# set LENGTH using an environment variable; the degree is ~ 287 * LENGTH 

sha1random::  
	$(clangs) -DPARAM=1 -DLENGTH=$(LENGTH) sha.c -o build/sha.3.s
	$(base) \
--load build/sha.3.s --compile sha1n $(LENGTH) --keygen \
--input Msg random --eval sha1n $(LENGTH) 

sha256random::  
	$(clangs) -DPARAM=1 -DLENGTH=$(LENGTH) sha.c -o build/sha.3.s
	$(base) \
--load build/sha.3.s --compile sha256n $(LENGTH) --keygen \
--input Msg random --eval sha256n $(LENGTH) 


# ---------------------------------------------------------------------------
# diverse tests in trivial.c (in CLI style)
trivials::
	(for p in 1 2 3 4 5 6 7; do make PARAM=$$p OPT=-q trivial ; done) | grep test

build/trivial.$(PARAM).exe : build pinocchio.c trivial.c 
	$(clangw) -DPARAM=$(PARAM) pinocchio.c trivial.c -o $@

trivial: build/trivial.$(PARAM).exe build/trivial.$(PARAM).s 
	@echo Testing trivial -DPARAM=$(PARAM)
	build/trivial.$(PARAM).exe
	$(base) --load build/trivial.$(PARAM).s --compile outsource 2 --keygen --input In 3,5,7 --eval outsource 2 --output Out

# testing a small part of the C semantics: structures, globals, comparisons, ...
unit-test:: unit-test.wp

# ---------------------------------------------------------------------------
# testing load/save, re-using saved keys and inputs

# in CLI-style, we always reload both keys, as they are both used by --eval.
Name=two-matrix.2
bankIn=$(Name).In.bank
bankOut=$(Name).Out.bank
bankFile:: build/$(Name).s
	$(base) --load build/$(Name).s --compile outsource --keygen --save both $(Name) 
	$(base) --resume both $(Name) --input Input random --output Input to $(bankIn) --eval outsource 2 --output Output
	$(base) --resume both $(Name) --input Input from $(bankIn)                     --eval outsource 2 --output Output
	ls -l build/$(bankIn)

# in whole-program style, we can separately run each stage; 
# note the schedule info is implicitly saved and restored 
stages:: build/matrix.s matrix.wp 
	@echo 
	@echo ----- compile and keygen
	$(base) -w --format Simple --load build/matrix.s --schedule build/matrix.trace --emit build/matrix --compile --keygen --save prover Prover --save verifier Verifier
	@echo 
	@echo ----- prove
	$(base) -w --format Rich --resume prover Prover --prove 
	@echo 
	@echo ----- verify
	$(base) -w --format Rich --resume verifier Verifier --verify


# ---------------------------------------------------------------------------
# a simple test for branching & energy-saving

branch:: build/branch.s
	$(base) --load build/branch.s --compile f 0 --compile f 1 --keygen \
--input In 0,4 --eval f 0 --output Out \
--input In 1,4 --eval f 0 --output Out \
--input In 0,4 --eval f 1 --output Out \
--input In 1,4 --eval f 1 --output Out
	$(base) --load build/branch.s --compile g --keygen --input In 0,1 --eval g --output Out


# ---------------------------------------------------------------------------
# basic tests for field operations and QAPed crypto

field:: build/field.s
	$(base) --load build/field.s --compile outsource --keygen --input Input random --eval outsource


# ---------------------------------------------------------------------------
# Whole-program proving, verifying, and bootstrapping 

.SUFFIXES: .wp

.c.wp:
	$(clangw) $(COPT) -DWHOLE -DNATIVE -DSCHEDULE $*.c -o build/$*.exe
	$(clangs) $(COPT) -DWHOLE -DINTERPRET $*.c -o build/$*.s
	@peflags -x0x4000000 ./build/$*.exe > /dev/null # increase stack size for large samples
	./build/$*.exe
	$(base) --encoding $(INNER) -w --format Simple --load build/$*.s --schedule build/$*.trace --emit build/$* --compile --keygen --prove --verify
# note that we require Simple format for bootstrapping

# note that users of pinocchio.h must explicitly load pinocchio.c for
# native compilation; see our two running examples:
# - make whole.wp whole.verify
# - make wide-loop.wp wide-loop.verify (using different sizes)

clangboot = $(clangs) -DVERIFY -DINTERPRET_BOOT -I $(QAP_ARITH) 
baseboot = $(base) -w --encoding $(OUTER) --bootstrap 
# no key saving yet, otherwise keys are overwritten!


## Compile native executable versions of the (bootstrappable) verifier
## using VS, to take advantage of Windows' intrinsics in arith;
## this native code is much slower than its boot-level verification.

build/whole.exe: whole.c geppetto.h build ../windows/basic-tests.vcxproj $(QAP_ARITH)/*h $(QAP_ARITH)/*c
	$(msbuild) /p:Configuration=$(BLT),Platform=x64,GeppettoMacro1=WHOLE,GeppettoMacro2=VERIFY,GeppettoMacro3=$(TINY) ../windows/basic-tests.vcxproj /t:"Clean;Rebuild"
	mv ../bin/$(BLT)/basic-tests.exe build/whole.exe

build/wide-loop.exe: wide-loop.c geppetto.h build ../windows/basic-tests.vcxproj $(QAP_ARITH)/*h $(QAP_ARITH)/*c
	$(msbuild) /p:Configuration=$(BLT),Platform=x64,GeppettoMacro1=WIDELOOP,GeppettoMacro2=VERIFY,GeppettoMacro3=$(TINY) ../windows/basic-tests.vcxproj /t:"Clean;Rebuild"
	mv ../bin/$(BLT)/basic-tests.exe build/wide-loop.exe

whole.verify: whole.c
	$(clangboot) whole.c -o build/whole-verify.s 
	$(baseboot) --load build/whole-verify.s --schedule build/whole-verify.trace --emit build/whole-verify --compile --keygen --prove --verify

whole-test: whole.wp whole.verify


wide-loop.verify: wide-loop.c
	$(clangboot) wide-loop.c -o build/wide-loop-verify.s 
	$(baseboot) --load build/wide-loop-verify.s --schedule build/wide-loop-verify.trace --emit build/wide-loop-verify --compile --keygen --prove --verify

wide-loop-test: wide-loop.wp wide-loop.verify
