(** Bzip2 interface *)

(** The module [Bz] provides a very basic interface to the [bzip2]
  compression library. *)

(** {2 Datatypes & exceptions} *)

type in_channel
type out_channel

exception IO_error of string
(** Exception [IO_error] is raised when there is an error reading or
  writing on a compressed channel ; the string argument is the message 
  reported by the OS. *)

exception Data_error
(** Exception [Data_error] is raised when a data integrity error
   is detected during decompression. *)

exception Unexpected_EOF
(** Exception [Unexpected_EOF] is raised when an [in_channel]
   finishes before the logical end of stream is detected. *)

(** When any of these exception is raised, the channel is automatically closed
   (but you still have to close the Pervasives channel). *)


val version : string

(** {2 Input functions} *)

(** [open_in] opens a compressed stream that is located on a 
   [Pervasives] channel. 
   @param small  defaults to [false] ; when [true] it uses a different 
   method for decompressing that is slower but uses less memory. 
*)
external open_in : ?small:bool -> ?unused:string -> Pervasives.in_channel -> in_channel 
  = "mlbz_readopen"

(** [read] reads [len] characters in a string buffer. 
   @return number of bytes actually read, (a value less than [len] 
   means end of stream). 
   @raise End_of_file if end of stream was already reached. *)
external read : in_channel -> buf:string -> pos:int -> len:int -> int
  = "mlbz_read"

(** If there's some data after the compressed stream that you want to read 
  from the same [Pervasives] [in_channel], 
   use [read_get_unused]. *)
external read_get_unused : in_channel -> string
  = "mlbz_readgetunused"

external close_in : in_channel -> unit
  = "mlbz_readclose"

(** {2 Output funcions} *)

(** [open_out] creates an [out_channel] from a [Pervasives]
   channel. Once the write operations are finished and the compressed channel 
   is closed, it is possible to continue writing on the [Pervasives]
   channel. However, reading back requires special care (cf. above).
   @param block block size to use for compresion. It is a value between 1 
   and 9 inclusive. 9 is the default and provides best compression but 
   takes most memory. *)
external open_out : ?block:int -> Pervasives.out_channel -> out_channel
  = "mlbz_writeopen"

external write : out_channel -> buf:string -> pos:int -> len:int -> unit
  = "mlbz_write"

external close_out : out_channel -> unit
  = "mlbz_writeclose"


(** {2 In-memory compression} *)

(** These functions compress & decompress to and from string buffers. *)

external compress : ?block:int -> string -> pos:int -> len:int -> string
  = "mlbz_compress"

external uncompress : ?small:bool -> string -> pos:int -> len:int -> string
  = "mlbz_uncompress"
