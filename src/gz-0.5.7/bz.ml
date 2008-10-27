type in_channel
type out_channel

exception IO_error of string
exception Data_error
exception Unexpected_EOF

let _ = begin
  Callback.register_exception "mlbz_io_exn" (IO_error "") ;
  Callback.register_exception "mlbz_data_exn" Data_error ;
  Callback.register_exception "mlbz_eof_exn" Unexpected_EOF
end

external library_version : unit -> string
  = "mlbz_version"

let version = library_version ()

external open_in : ?small:bool -> ?unused:string -> Pervasives.in_channel -> in_channel 
  = "mlbz_readopen"

external read : in_channel -> buf:string -> pos:int -> len:int -> int
  = "mlbz_read"

external read_get_unused : in_channel -> string
  = "mlbz_readgetunused"

external close_in : in_channel -> unit
  = "mlbz_readclose"

external open_out : ?block:int -> Pervasives.out_channel -> out_channel
  = "mlbz_writeopen"

external write : out_channel -> buf:string -> pos:int -> len:int -> unit
  = "mlbz_write"

external close_out : out_channel -> unit
  = "mlbz_writeclose"


external compress : ?block:int -> string -> pos:int -> len:int -> string
  = "mlbz_compress"

external uncompress : ?small:bool -> string -> pos:int -> len:int -> string
  = "mlbz_uncompress"
