/**
@file

Copyright John Reid 2007, 2008, 2012

*/

#include "hmm-pch.h"

#include <boost/python.hpp>
#include "hmm/defs.h"

#include <myrrh/fpu.h>

#include <float.h>

int major_version() { return 1; }
int minor_version() { return 0; }
std::string version() { return HMM_MAKE_STRING( major_version() << "." << minor_version() ); }

unsigned denormal_control_mask() { return _MCW_DN; }
unsigned denormal_save() { return _DN_SAVE; }
unsigned denormal_flush() { return _DN_FLUSH; }
	
unsigned interrupt_exception_mask() { return _MCW_EM; }
unsigned interrupt_exception_invalid() { return _EM_INVALID; }
unsigned interrupt_exception_denormal() { return _EM_DENORMAL; }
unsigned interrupt_exception_zero_divide() { return _EM_ZERODIVIDE; }
unsigned interrupt_exception_overflow() { return _EM_OVERFLOW; }
unsigned interrupt_exception_underflow() { return _EM_UNDERFLOW; }
unsigned interrupt_exception_inexact() { return _EM_INEXACT; }

unsigned infinity_control_mask() { return _MCW_IC; }
unsigned infinity_control_affine() { return _IC_AFFINE; }
unsigned infinity_control_projective() { return _IC_PROJECTIVE; }

unsigned rounding_control_mask() { return _MCW_RC; }
unsigned rounding_control_chop() { return _RC_CHOP; }
unsigned rounding_control_up() { return _RC_UP; }
unsigned rounding_control_down() { return _RC_DOWN; }
unsigned rounding_control_near() { return _RC_NEAR; }
	
unsigned precision_control_mask() { return _MCW_PC; }
unsigned precision_control_24() { return _PC_24; }
unsigned precision_control_53() { return _PC_53; }
unsigned precision_control_64() { return _PC_64; }

///** Returns old state. see http://msdn2.microsoft.com/en-us/library/e9b52ceh(VS.80).aspx */
//int control87_2( unsigned new_state, unsigned mask, bool set_x86_cw, unsigned x86_cw, bool set_sse2_cw, unsigned sse2_cw )
//{
//	return _control87_2( new_state, mask, set_x86_cw ? &x86_cw : 0, set_sse2_cw ? &sse2_cw : 0 );
//}


BOOST_PYTHON_MODULE( _fpu )
{
	using namespace boost::python;
	using boost::python::arg;

	scope().attr("__doc__") = 
		"Gets/sets various bits of the floating point control word\r\n"
		"See U{http://msdn2.microsoft.com/en-us/library/e9b52ceh(VS.80).aspx}";

	def( 
		"__version__",
		version,
		"The version of the C++ fpu python extension" );

	def( 
		"summarise_status",
		myrrh::fpu_control_word_summary,
		"A string stating the FPU status." );

	def( 
		"control_fp",
		_controlfp,
		"Sets the fpu control state. Returns the old state.",
		( arg("new_state")=0, arg("mask")=0 ) );

	def( 
		"control_87",
		_control87,
		"Sets the fpu control state. Returns the old state.",
		( arg("new_state")=0, arg("mask")=0 ) );

//	def( 
//		"control_87_2",
//		control_87_2,
//		"Sets the fpu control state. Returns the old state.",
//		( arg("new_state")=0, arg("mask")=0, arg("set_x86_cw")=false, arg("x86_cw")=0, arg("set_sse2_cw")=false, arg("sse2_cw")=0 ) );

	def( "denormal_control_mask", denormal_control_mask );
	def( "denormal_save", denormal_save );
	def( "denormal_flush", denormal_flush );

	def( "interrupt_exception_mask", interrupt_exception_mask );
	def( "interrupt_exception_invalid", interrupt_exception_invalid );
	def( "interrupt_exception_denormal", interrupt_exception_denormal );
	def( "interrupt_exception_zero_divide", interrupt_exception_zero_divide );
	def( "interrupt_exception_overflow", interrupt_exception_overflow );
	def( "interrupt_exception_underflow", interrupt_exception_underflow );
	def( "interrupt_exception_inexact", interrupt_exception_inexact );

	def( "infinity_control_mask", infinity_control_mask );
	def( "infinity_control_affine", infinity_control_affine );
	def( "infinity_control_projective", infinity_control_projective );

	def( "rounding_control_mask", rounding_control_mask );
	def( "rounding_control_chop", rounding_control_chop );
	def( "rounding_control_up", rounding_control_up );
	def( "rounding_control_down", rounding_control_down );
	def( "rounding_control_near", rounding_control_near );

	def( "precision_control_mask", precision_control_mask );
	def( "precision_control_24", precision_control_24 );
	def( "precision_control_53", precision_control_53 );
	def( "precision_control_64", precision_control_64 );
}


