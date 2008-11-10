type color = int
let black = 1
let current_file = ref None
let contents = ref []
let resources = ref []
let max_x = ref 0.
let max_y = ref 0.
let current_x = ref 0.
let current_y = ref 0.

let set_max () =
    max_x := max !max_x !current_x;
    max_y := max !max_y !current_y

let clear_contents () =
    max_x := 0.;
    max_y := 0.;
    current_x := 0.;
    current_y := 0.;
    contents := [];
    resources := []

let pages = ref [] 

let open_file (x : string) = 
    clear_contents ();
    pages := [];
    current_file := Some x

let open_graph size =
    clear_contents ()

let add_page () = 
    let page = 
        let font =
            Pdf.Dictionary
            [("/Type", Pdf.Name "/Font");
            ("/Subtype", Pdf.Name "/Type1");
            ("/BaseFont", Pdf.Name "/Times-Italic")]
        in
        let contents = !contents in
        {(Pdfdoc.blankpage (Units.PdfPoint, !max_x, !max_y)) with
            Pdfdoc.content = contents;
            resources = Pdf.Dictionary (("/Font", Pdf.Dictionary [("/F0",
            font)]) :: !resources)}
    in
    pages := page :: !pages;
    clear_contents ()

let close_graph () =
    add_page ();
    let pdf = 
        let empty = 
            { Pdf.empty with Pdf.minor = 5 } 
        in
            let pdf, pageroot =
                Pdfdoc.add_pagetree (List.rev !pages) empty
            in
            Pdfdoc.add_root pageroot [] pdf
    in
    match !current_file with
    | Some "" 
    | None -> failwith "No file opened"
    | Some x ->
            Pdfwrite.pdf_to_file pdf x;
            clear_contents ()


let add_ops ops =
    let ops = (Pdfpages.stream_of_ops ops) in
    contents :=  ops :: !contents

let text_size string = 
    ((float_of_int (String.length string)) *. 5.), 6.

let last_x_string = ref 0.
let last_y_string = ref 0.

let draw_string string =
    if String.length string = 0 then ()
    else
        let x, y = text_size string in
        let ops = 
            [ Pdfpages.Op_cm (Transform.matrix_of_transform [Transform.Translate
            (!current_x, !current_y)]);
                Pdfpages.Op_BT; 
            Pdfpages.Op_Tf ("/F0", 12.);
            Pdfpages.Op_Tj string;
            Pdfpages.Op_ET;
            Pdfpages.Op_cm (Transform.matrix_of_transform [Transform.Translate
                        ((-. !current_x), (-. !current_y))]);

            ]
        in
        add_ops ops;
        current_x := !current_x +. (x *. 2.);
        current_y := !current_y +. (y *. 4.);
        set_max ()

let foreground = black

let lineto x y =
    let x = float_of_int x
    and y = float_of_int y in
    let ops = 
        [ 
            Pdfpages.Op_m (!current_x, !current_y);
            Pdfpages.Op_l (x, y);
            Pdfpages.Op_h;
            Pdfpages.Op_S]
    in
    add_ops ops;
    current_x := x;
    current_y := y; 
    set_max ()

let polyline lst =
    let ops = 
        match lst with
        | (a, b) :: tl ->
                (Pdfpages.Op_m (float_of_int a, float_of_int b)) ::
                    (List.fold_right (fun (a, b) acc ->
                        let a, b = float_of_int a, float_of_int b in
                        current_x := a;
                        current_y := b;
                        (Pdfpages.Op_l (a, b)) :: acc) tl [Pdfpages.Op_S])
        | _ -> failwith "Polyline needs at least two points"
    in
    add_ops ops;
    set_max ()

let move_to x y =
    current_x := x;
    current_y := y;
    set_max ()

let plot x y = 
    move_to (float_of_int x) (float_of_int y)

let red = black
let set_color _ = ()
let size_x () = int_of_float !max_x
let size_y () = int_of_float !max_y

let display () = ()

let text_size string = 
    let a, b = text_size string in
    int_of_float a, int_of_float b

let moveto a b =
    move_to (float_of_int a) (float_of_int b)
