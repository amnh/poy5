#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <caml/mlvalues.h>
#include <caml/fail.h>
#include <caml/alloc.h>
#include <caml/callback.h>

#include <bzlib.h>
#include "io.h"

#ifdef BZ_PREFIX
#define BZ_P(a) BZ2_ ## a
#else
#define BZ_P(a) a
#endif 


value mlbz_readclose(value chan) ;
value mlbz_writeclose(value chan) ;

/* ERROR REPORTING */
/* checks bzerror and raises the corresponding Caml exception */
/* then closes the channel */
static void mlbz_error(int bzerror, char *msg, value chan, int read_str)
{
  if(bzerror<0){
    if(read_str)
      mlbz_readclose(chan) ;
    else
      mlbz_writeclose(chan) ;
  }
  switch(bzerror){
  case BZ_PARAM_ERROR :
    invalid_argument(msg) ; break ;
  case BZ_MEM_ERROR : 
    raise_out_of_memory() ; break ;
  case BZ_DATA_ERROR :
  case BZ_DATA_ERROR_MAGIC :
    raise_constant(*caml_named_value("mlbz_data_exn")) ; break ;
  case BZ_IO_ERROR :
    raise_with_string(*caml_named_value("mlbz_io_exn"), 
		      strerror(errno)) ; break ;
  case BZ_UNEXPECTED_EOF :
    raise_constant(*caml_named_value("mlbz_eof_exn")) ; break ;
  default :
    failwith(msg) ;
  }
}

/* MEMORY HANDLING */ 
/* The Caml abstract_value contains 3 values : */
/* - BZFILE* pointer (compressed stream) */
/* - FILE* pointer (underlying C stream) */
/* - bool value (flag indicating an end of stream */
static value Val_twoptr(FILE * ptr1, BZFILE * ptr2)
{
  value result;
  result= alloc_small(3, Abstract_tag);
  Field(result, 0) = (value) ptr1 ; 
  Field(result, 1) = (value) ptr2 ; 
  Field(result, 2) = Val_false ;     /* flag for the eof state */
  return result;
}

#define Cfile_val(v)  ((FILE*)(Field((v), 0)))
#define Bzfile_val(v) ((BZFILE*)(Field((v), 1)))
#define Bz_eof(v)     (Bool_val(Field((v), 2)))

static void set_eof_flag(value chan)
{
  Field(chan, 2) = Val_true ;
}

/* converts a Caml channel to a C FILE* stream */
static FILE * stream_of_channel(value chan, const char * mode)
{
  int des ;
  FILE * res ;
  struct channel *c_chan = Channel(chan) ;
  if(c_chan==NULL)
    return NULL;
  des = dup(c_chan->fd) ;
  res = fdopen(des, mode) ;
  return res ;
}

value mlbz_version(value unit)
{
  return copy_string(BZ_P(bzlibVersion())) ;
}

/* INPUT FUNCTIONS */
value mlbz_readopen(value small, value unused, value chan)
{
  int bzerror ;
  int c_small = 0 ;
  char *c_unused = NULL ;
  int c_nunused = 0;
  FILE * c_chan ;
  BZFILE * bz_chan ;

  if(Is_block(small))
    c_small = Bool_val(Field(small, 0)) ;
  if(Is_block(unused)){
    c_unused  = String_val(Field(unused, 0)) ;
    c_nunused = string_length(Field(unused, 0));
  }
  c_chan = stream_of_channel(chan, "rb") ;
  bz_chan = BZ_P(bzReadOpen)(&bzerror, c_chan, c_small, 0, c_unused, c_nunused) ;
  mlbz_error(bzerror, "Bz.open_in", chan, 1) ;
  return Val_twoptr(c_chan, bz_chan) ;
}

value mlbz_read(value chan, value buf, value pos, value len)
{
  int res ;
  int bzerror ;
  char *c_buf ;
  int c_pos = Int_val(pos);
  int c_len = Int_val(len);
  
  if(Bz_eof(chan))
    raise_end_of_file() ;    
  if((c_len + c_pos > string_length(buf))
     || (c_len < 0) || (c_pos < 0))
    invalid_argument("Bz.read") ;
  c_buf = String_val(buf) + c_pos ;
  res = BZ_P(bzRead)(&bzerror, Bzfile_val(chan), c_buf, c_len) ;
  if(bzerror == BZ_STREAM_END)
    set_eof_flag(chan) ;
  mlbz_error(bzerror, "Bz.read", chan, 1) ;
  return Val_int(res) ;
}

value mlbz_readclose(value chan)
{
  int bzerror ;
  BZ_P(bzReadClose)(&bzerror, Bzfile_val(chan)) ;
  fclose(Cfile_val(chan)) ;
  {
    void **p = (void **) chan;
    p[0] = NULL;
    p[1] = NULL;
  }
  /*  mlbz_error(bzerror, "Bz.close_in", chan, 1); */
  return Val_unit ;
}

value mlbz_readgetunused(value chan)
{
  value result;
  int bzerror ;
  void *unused ;
  int nunused ;
  BZ_P(bzReadGetUnused)(&bzerror, Bzfile_val(chan), &unused, &nunused) ;
  mlbz_error(bzerror, "Bz.read_get_unused: not at end of stream", chan, 1) ;
  result = alloc_string(nunused) ;
  memcpy(String_val(result), unused, nunused) ;
  return result;
}


/* OUTPUT FUNCTIONS */
value mlbz_writeopen(value block, value chan)
{
  int bzerror ;
  int c_block = 9 ;
  FILE * c_chan ;
  BZFILE * bz_chan ;
  if(Is_block(block))
    c_block = Int_val(Field(block, 0)) ;
  c_chan = stream_of_channel(chan, "wb");
  bz_chan = BZ_P(bzWriteOpen)(&bzerror, c_chan, c_block,
			0, 0) ;
  mlbz_error(bzerror, "Bz.open_out", chan, 0);
  return Val_twoptr(c_chan, bz_chan);
}

value mlbz_write(value chan, value buf, value pos, value len)
{
  int bzerror ;
  char *c_buf ;
  int c_pos = Int_val(pos) ;
  int c_len = Int_val(len) ;
  if((c_len + c_pos > string_length(buf))
     || (c_len < 0) || (c_pos < 0))
    invalid_argument("Bz.write") ;
  c_buf = String_val(buf) + c_pos ;
  BZ_P(bzWrite)(&bzerror, Bzfile_val(chan), c_buf, c_len) ;
  mlbz_error(bzerror, "Bz.write", chan, 0) ;
  return Val_unit;
}

value mlbz_writeclose(value chan)
{
  int bzerror ;
  BZ_P(bzWriteClose)(&bzerror, Bzfile_val(chan), 0, NULL, NULL) ;
  fclose(Cfile_val(chan)) ;
  {
    void **p = (void **) chan;
    p[0] = NULL;
    p[1] = NULL;
  }
  /*  mlbz_error(bzerror, "Bz.close_out", chan, 0); */
  return Val_unit ;
}


/* IN MEMORY COMPRESSION */
value mlbz_compress(value block, value src, value pos, value len)
{
  value result ;
  int c_block=9 ;
  int c_pos = Int_val(pos) ;
  int c_len = Int_val(len) ;
  int r ;
  unsigned int dst_buf_len, dst_len ;
  char *src_buf, *dst_buf;

  if(Is_block(block))
    c_block = Int_val(Field(block, 0)) ;
  if(c_block < 1 || c_block > 9
     || c_pos < 0 || c_len < 0
     || c_pos + c_len > string_length(src))
    invalid_argument("Bz.compress") ;
  src_buf = String_val(src) + c_pos;
  dst_buf_len = c_len * 1.01 + 600 ;
  dst_buf = malloc(dst_buf_len) ;
  if(dst_buf == NULL)
    raise_out_of_memory();
  while(1) {
    dst_len = dst_buf_len;
    r = BZ_P(bzBuffToBuffCompress)(dst_buf, &dst_len, src_buf, c_len, c_block, 0, 0) ;
    if(r == BZ_OK) {
      break;
    } else if(r == BZ_OUTBUFF_FULL){
      char *new_buf;
      dst_buf_len *= 2;
      new_buf = realloc(dst_buf, dst_buf_len);
      if(new_buf == NULL) {
	free(dst_buf);
	raise_out_of_memory();
      }
      dst_buf = new_buf;
    } else {
      free(dst_buf);
      raise_out_of_memory();
    }
  }
  result = alloc_string(dst_len);
  memcpy(String_val(result), dst_buf, dst_len);
  free(dst_buf);
  return result ;
}

value mlbz_uncompress(value small, value src, value pos, value len)
{
  value result ;
  int c_small = 0 ;
  int c_pos = Int_val(pos) ;
  int c_len = Int_val(len) ;
  int r ;
  unsigned int dst_buf_len, dst_len ;
  char *src_buf, *dst_buf;

  if(Is_block(small))
    c_small = Bool_val(Field(small, 0)) ;
  if(c_pos < 0 || c_len < 0
     || c_pos + c_len > string_length(src))
    invalid_argument("Bz.uncompress") ;
  src_buf = String_val(src) + c_pos;
  dst_buf_len = c_len * 2 ;
  dst_buf = malloc(dst_buf_len) ;
  if(dst_buf == NULL)
    raise_out_of_memory();
  while(1) {
    dst_len = dst_buf_len;
    r = BZ_P(bzBuffToBuffDecompress)(dst_buf, &dst_len, src_buf, c_len, c_small, 0) ;
    if(r == BZ_OK) 
      break;
    else 
      switch(r) {
      case BZ_OUTBUFF_FULL : {
	char *new_buf;
	dst_buf_len *= 2;
	new_buf = realloc(dst_buf, dst_buf_len);
	if(new_buf == NULL) {
	  free(dst_buf);
	  raise_out_of_memory();
	}
	dst_buf = new_buf;
      } break ;
      case BZ_MEM_ERROR :
	free(dst_buf) ;
	raise_out_of_memory() ;
      case BZ_DATA_ERROR :
      case BZ_DATA_ERROR_MAGIC :
	raise_constant(*caml_named_value("mlbz_data_exn")) ; 
      case BZ_UNEXPECTED_EOF :
	raise_constant(*caml_named_value("mlbz_eof_exn")) ;
      }
  }
  result = alloc_string(dst_len);
  memcpy(String_val(result), dst_buf, dst_len);
  free(dst_buf);
  return result ;
}
