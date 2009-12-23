module Parsed : sig
    val of_channel : in_channel -> string -> Nexus.File.nexus
end

module File : sig
    type characters = 
        | Range of (int * int)
        | Single of int
        | All

    type cost_change = (bool * int * (string * string) * int)

    type char_change =
        | Additive of characters list
        | NonAdditive of characters list
        | Active of characters list
        | Inactive of characters list
        | Sankoff of characters list
        | Weight of (int * characters list)

    type char_name = string
    type gappy = [ `Nogap | `Gap ] option
    type command = 
        | Nstates of [ `Dna of gappy | `Protein of gappy | `Number of int ] option
        | Ccode of char_change list
        | Cost of cost_change list
        | Tread of string
        | Xread of string
        | Ignore
        | Charname of char_name list


    val is_hennig : FileStream.f -> bool
end

module Grammar : sig
    type token
    val xread : (Lexing.lexbuf -> token) -> Lexing.lexbuf -> (int * int * string)
    val command : (Lexing.lexbuf -> token) -> Lexing.lexbuf -> File.command
end

module Lexer : sig
    exception Eof
    val token : Lexing.lexbuf -> Grammar.token
    val xread : Lexing.lexbuf -> Grammar.token
end
