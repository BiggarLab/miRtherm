CXX ?= g++
CC ?= gcc
CFLAGS ?= -O2
VSRC = core/vienna/cofold.c core/vienna/params.c core/vienna/energy_par.c core/vienna/fold_vars.c core/vienna/utils.c core/vienna/globals.c
VOBJ = $(VSRC:.c=.o)

mirtherm: cli/main.cpp core/mirtherm.cpp core/mirtherm.h $(VOBJ)
	$(CXX) $(CFLAGS) -std=c++11 -o $@ cli/main.cpp core/mirtherm.cpp $(VOBJ) -lm

core/vienna/%.o: core/vienna/%.c
	$(CC) $(CFLAGS) -std=gnu89 -fgnu89-inline -w -Icore/vienna -c $< -o $@

clean:
	rm -f mirtherm $(VOBJ)

# WebAssembly build for the dashboard (needs Emscripten)
wasm: docs/mirtherm_wasm.js
docs/mirtherm_wasm.js: wasm/wrapper.cpp core/mirtherm.cpp core/mirtherm.h $(VSRC)
	emcc -O2 -std=gnu89 -fgnu89-inline -w -Icore/vienna -c $(VSRC)
	em++ -O2 -std=c++11 -fno-exceptions -o $@ wasm/wrapper.cpp core/mirtherm.cpp \
	  cofold.o params.o energy_par.o fold_vars.o utils.o globals.o \
	  -s MODULARIZE=1 -s EXPORT_NAME=MirthermModule -s ALLOW_MEMORY_GROWTH=1 \
	  -s ENVIRONMENT=web,worker,node \
	  -s EXPORTED_FUNCTIONS='["_mt_configure","_mt_header","_mt_scan","_malloc","_free"]' \
	  -s EXPORTED_RUNTIME_METHODS='["cwrap","UTF8ToString"]'
	rm -f cofold.o params.o energy_par.o fold_vars.o utils.o globals.o
