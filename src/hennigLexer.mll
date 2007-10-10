(* Hennig lexer *)
{
    open HennigParser
    exception Eof

    let ( --> ) a b = b a

    let keyword_table = Hashtbl.create 53

    let token_table = 
        let make_all_prefixes min_len fullcommand output acc =
            let max_len = String.length fullcommand in 
            let rec prepend len acc =
                if len > max_len then acc
                else 
                    prepend (len + 1) (((String.sub fullcommand 0 len), output) :: acc)
            in
            prepend min_len acc
        in
        [] 
        --> make_all_prefixes 2 "ccode" CCODE
        --> make_all_prefixes 2 "costs" COST
        --> make_all_prefixes 2 "proc" PROCESS
        --> make_all_prefixes 2 "optcode" OPTCODE
        --> make_all_prefixes 2 "cnames"    CHARNAMECMD
        --> make_all_prefixes 2 "nstates" NSTATES
        --> (fun x -> ("dna", DNA) :: ("prot", PROTEINS) :: ("num", NUMBER) ::
            x)

    let _ = 
        List.iter (fun (keyw, tok) -> 
            Hashtbl.add keyword_table (String.uppercase keyw) tok) 
        token_table

    let is_prefix a b =
        let la = String.length a in
        la <= (String.length b) &&
        (a = (String.sub b 0 la))
}

rule token = parse
      [ ' ' '\t' '\n' '\010' '\013' '\012' ]    { token lexbuf }
    | [';'] { SEMICOLON }
    | [ '-' ] { DASH }
    | [ '0' - '9']+ as id { INT id }
    | [ '[' ] { LSQ }
    | [ ']' ] { RSQ }
    | [ '{' ] { characters_names lexbuf }
    | [ '+' ] { PLUS }
    | [ '*' ] { STAR }
    | [ '/' ] { BACKSLASH }
    | [ '.' ] { DOT }
    | [ '>' ] { GT }
    | [ '=' ] { EQUAL }
    | [ '?' ] { QUESTION }
    | [ '(' ] { LPARENT }
    | [ ')' ] { RPARENT }
    | [ 'a'-'z' 'A'-'Z']+ as word { 
        let uword = String.uppercase word in
        try Hashtbl.find keyword_table uword with 
        | Not_found -> 
                if is_prefix uword "TREAD" then rawtree lexbuf
                else if is_prefix uword "XREAD" then raw lexbuf
                else WORD word }
    | [ ^ '\000' ] as ch { CHAR ch } 
    | eof       { raise Eof }
and characters_names = parse
     [^ ';']+[';'] as r { CHARNAME r }
and raw = parse
     [^ ';']* as d      { DATA d }
and rawtree = parse
      [ ' ' '\t' '\n' '\010' '\013' '\012' ]    { rawtree lexbuf }
    | [ '\'' ]  { ignore_quote2 lexbuf }
    | [^ '\'' ';']* as d      { TREES d }
and ignore_quote2 = parse
    | [^ '\'']+['\'']   { rawtree lexbuf }
and ignore_quote = parse
    | [^ '\'']+['\'']   { xread lexbuf }
and xread = parse 
      [ ' ' '\t' '\n' '\010' '\013' '\012' ]    { xread lexbuf }
    | [ '\'' ]  { ignore_quote lexbuf }
    | [ '0'-'9']+ as integer { INT integer }
    | ['a'-'z' 'A'-'Z' '_' '-' '.' '?'][^ '\000']+ eof as data    { DATA data }
