#ifndef _stdfunc_h
#define _stdfunc_h

#define SAFE_DELETE(x) {if (x) {delete (x); (x)=NULL;}}
#define SAFE_DELETE_ARRAY(x) {if (x) {delete [] (x); (x)=NULL;}}

#endif
