open Printf
let printf format = printf format
let failwithf format = ksprintf (failwith) format


let gap = 4 
let gapcost = 2. 
let cost_matrix = [|
    [|  0.; 1.;  1.; 1.; gapcost|];
    [|  1.; 0.;  1.; 1.; gapcost|];
    [|  1.; 1.;  0.; 1.; gapcost|];
    [|  1.; 1.;  1.; 0.; gapcost|];
    [|  gapcost; gapcost; gapcost; gapcost; 0.|];
|]

type dir      = Root | Align | Delete | Insert


let choose_dir dirs =
             if List.mem Delete dirs then Delete
        else if List.mem Align  dirs then Align
        else if List.mem Insert dirs then Insert
        else match dirs with | [Root] -> Root | _ -> assert false

let mem = ref None

let create_mem a b =
        let a = max a b in
        let res = Array.make_matrix a a (~-.1.0,(0,[Root])) in
        mem := Some res;
        res

let get_mem a b =    
    create_mem ((Array.length a)) ((Array.length b)) 

let dirToString = function
        | Root   -> "X" | Align  -> "\\"
        | Delete -> "|" | Insert -> "-"
    
let print_mem mem oc =
        for i = 0 to (Array.length mem)-1 do
        for j = 0 to (Array.length mem.(i))-1 do
            let cost,(indel,dirs) = mem.(i).(j) in
            fprintf oc "|(%d,%d) %s %02d[%02d] "
                i j (dirToString (choose_dir dirs)) (int_of_float cost) indel;
        done;
        fprintf oc "|\n%!";
    done;
    fprintf oc "|\n%!"
    
(* minimum of three with annotations; this will work for all methods *)
let min3_ a at b bt c ct =
    if a > b then begin
        if b > c then (c,[ct])
        else if b = c then ((min b c),[bt;ct])
        else (b,[bt])
    end else if a = b then begin
        if b = c then ((min a (min b c)),[at;bt;ct])
        else if a > c then (c,[ct])
        else ((min a b),[at;bt])
    end else begin
        if c > a then (a,[at])
        else if c = a then ((min a c),[at;ct])
        else (c,[ct])
    end

(* [align_2 mem x y m t] Align the sequence [x] and [y] with the cost
     * matrix from [m] and branch length [t], and completely fills the matrix of
     * [mem]. Return minimum cost; the alignment can be recovered from [mem] *)
    let align_2 mem x y (* m:dyn_model *) (* tx:float) (ty:float *) oc =
        (*let xlen = Sequence.length x and ylen = Sequence.length y in
        Numerical.FPInfix.set_ops (2 * (xlen + ylen));
        let cost = let cf = cost m tx ty in (fun a b -> fst (cf a b) ) in
        let gap = Alphabet.get_gap m.alph in
        let get a b = Sequence.get a b in
        *)
        let xlen = Array.length x and ylen = Array.length y in
        assert( xlen <= ylen );
        fprintf oc "align 2 , lenX=%d, lenY=%d\n%!" xlen ylen;
        let cost addx addy = cost_matrix.(addx).(addy) in
        let get_cost i j =
            if i = 0 && j = 0 then begin
                (0.0,[Root])
            end else if i = 0 then begin
                (fst mem.(i).(j-1) +. (cost gap y.(j))),[Insert]
            end else if j = 0 then begin
                (fst (mem.(i-1).(j))) +. (cost x.(i) gap),[Delete]
            end else begin
                min3_ ((fst mem.(i-1).(j))   +. (cost x.(i) gap)) Delete
                      ((fst mem.(i).(j-1))   +. (cost gap y.(j))) Insert
                      ((fst mem.(i-1).(j-1)) +. (cost x.(i) y.(j))) Align
            end
        in
        for i = 0 to xlen - 1 do
            for j = 0 to ylen - 1 do
                let c,ds = get_cost i j in
                mem.(i).(j) <- c,(0,ds)
            done;
        done;
        (*Numerical.FPInfix.reset ();
        if debug_aln then begin
            let ex,ey = alignments mem x y m in
            printf "C1: "; print_s ex m.alph;
            printf "C2: "; print_s ey m.alph;
            printf "%f ME: " (fst mem.(xlen-1).(ylen-1));
            print_s (backtrace ~filter_gap:false m mem x y) m.alph;
        end;*)
        fst mem.(xlen-1).(ylen-1)


let () =
    let oc =  open_out "algnfull.out" in
    let filename = "ukkfile.in" in
    let lines = ref [] in
    let chan = open_in filename in
    try 
    while true do
    lines := input_line chan :: !lines
    done;
    with End_of_file ->
    close_in chan;
    let lines = List.rev !lines in
    fprintf oc "read in str=\n%!";
    List.iter (fun str -> fprintf oc "%s\n%!" str) lines;
    let seqlst = List.fold_right (fun seqstr acc ->
        if (String.length seqstr)>0 && 
        (seqstr.[0]='a'||seqstr.[0]='c'||seqstr.[0]='t'||seqstr.[0]='g') then
        begin
        let seqlst = ref [] in
        for i=0 to (String.length seqstr)-1 do
            seqlst := 
            (match seqstr.[i] with
            | 'a' -> 0::!seqlst
            | 'c' -> 1::!seqlst
            | 'g' -> 2::!seqlst
            | 't' -> 3::!seqlst
            | _ -> !seqlst )
            ;
        done;
        (Array.of_list (List.rev !seqlst))::acc
        end
        else acc
    ) lines [[||]] 
    in
    fprintf oc "seqarr =\n%!";
     List.iter (
        fun seqlst ->
            fprintf oc "[%!";
            Array.iter (fun x ->fprintf oc "%d" x) seqlst; 
        fprintf oc "]\n%!";
        ) seqlst; 
    let seqarr1,seqarr2 = List.hd seqlst, List.nth seqlst 1 in
    let seqarr1,seqarr2 = 
        Array.append [|gap|] seqarr1, 
        Array.append [|gap|] seqarr2
    in
    let seqarr1,seqarr2 = 
        if (Array.length seqarr1)<=(Array.length seqarr2) then
            seqarr1,seqarr2
        else 
            seqarr2,seqarr1
    in
    let ukkmem = get_mem seqarr1 seqarr2 in
    let cost = align_2 ukkmem seqarr1 seqarr2 oc in
    fprintf oc "float cost = %f\n%!" cost;
    printf "%d%!" (int_of_float cost); 
    close_out oc
    
