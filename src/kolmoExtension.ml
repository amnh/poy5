
let () =
    let module M = 
        Camlp4.Register.OCamlSyntaxExtension (KolmoGrammar.IdSKOcaml)
        (KolmoGrammar.SKOcamlLanguageExt) 
    in 
    ()
