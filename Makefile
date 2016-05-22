PROG_LIST = cytoscore configure_cytoscore
TEST_LIST = motor_test diffusion_test msd_analysis

.PHONY : default
default :
	@echo
	@echo "Usage: make <prog>"
	@echo "   where prog is one of: ($(PROG_LIST))"
	@echo

$(PROG_LIST): force-build
	cd src; $(MAKE) ../$@

$(TEST_LIST): force-build
	cd test; $(MAKE) $@

.PHONY : test
test:
	./test.sh

.PHONY : clean
clean :
	rm -f src/*.o test/*.o

.PHONY : clean-output
clean-output :
	rm -f *.posit *.config *.thermo *.initial_config.* *.final_config *.final_config.* \
	*.checkpoint.* *.checkpoint *.thermo_ext sphero.crosslinks.* test-log

.PHONY : clean-all
clean-all : clean clean-output
	rm -f $(PROG_LIST)

force-build:

