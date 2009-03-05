let cost = ref [||]
let vertical : float array array ref = ref [||]
let horizontal : float array array ref = ref [||]
let diagonal : float array array ref = ref [||]
let direction = ref [||]
let event = ref [||]
let size = ref 0

type matrix = {
    event_cost : float;
    first_event_cost : float;
    matrix : Cost_matrix.Two_D.m;
}

let minimum_event matrix apos bpos =
    (!event).(bpos).(apos) <- 
        matrix.event_cost +. (min 
            (min (!event).(bpos - 1).(apos - 1) (!event).(bpos -
            1).(apos))
                (!event).(bpos).(apos - 1))

let eventf event_cost y x matrix = 
    if event_cost <> 0. then matrix.first_event_cost
    else (!event).(x).(y) 

let eventf_ac event_cost y x matrix =
    if event_cost = 0. then 0.
    else (!event).(x).(y)

let costf matrix y x = (!cost).(x).(y)
let verticalf matrix y x = (!vertical).(x).(y)
let horizontalf matrix y x = (!horizontal).(x).(y)
let diagonalf matrix y x = (!diagonal).(x).(y)

let large_float = 100000.

let is_in a b = 0 <> (a land b)

let diag gap a b =
    if (is_in gap (a land b)) then 0.
    else large_float

let go gap a i opening =
    if i = 1 && (is_in gap (Sequence.get a i)) then 0.
    else if i > 1 && (is_in gap (Sequence.get a i)) && (not (is_in gap
    (Sequence.get a (i - 1)))) then 0.
    else opening

let subst matrix a b =
    float_of_int (Cost_matrix.Two_D.cost a b matrix.matrix)

let go' matrix opening gap a b = 
    (subst matrix a b) +.  (if not (is_in gap a) then 0. else opening)

let ge a gap extension = 
    if is_in gap a then 0. else float_of_int extension

let min4 a b c d =
    min (min a b) (min c d)

let select_best_direction apos bpos =
    let a = 1
    and h = 2
    and v = 4
    and d = 8 in
    let algn = (!cost).(bpos).(apos)
    and horz = (!horizontal).(bpos).(apos)
    and vert = (!vertical).(bpos).(apos)
    and diag = (!diagonal).(bpos).(apos) in
    (!direction).(bpos).(apos) <-
        if algn < horz then
            if algn < vert then 
                if algn < diag then a
                else if diag < algn then d
                else a lor d
            else if vert < algn then 
                if vert < diag then v
                else if diag < vert then d
                else v lor d
            else 
                if algn < diag then a
                else if diag < algn then a
                else a lor d
        else if horz < algn then 
            if horz < vert then
                if horz < diag then h
                else if diag < horz then d
                else h lor d
            else  if vert < horz then
                if vert < diag then v
                else if diag < vert then d
                else d lor v
            else 
                if vert < diag then v lor h
                else if diag < vert then d
                else d lor v lor h
        else 
            (* horz = algn *) 
            if horz < vert then
                if horz < diag then h lor a
                else if diag < horz then d
                else h lor a lor d
            else if vert < horz then 
                if vert < diag then v
                else if diag < vert then d
                else v lor d
            else 
                (* horz = vert *)
                if horz < diag then h lor v lor a
                else if diag < horz then d
                else 
                    (* horz = diag *)
                    h lor v lor a lor d

let affine_align gap opening matrix aseq apos bseq bpos =
    let abase = Sequence.get aseq apos
    and bbase = Sequence.get bseq bpos in
    let bextension = Cost_matrix.Two_D.cost bbase gap matrix.matrix
    and aextension = Cost_matrix.Two_D.cost abase gap matrix.matrix in
    (!cost).(bpos).(apos) <-
        min4
            ((!cost).(bpos - 1).(apos - 1) +. (subst matrix abase bbase))
            ((!diagonal).(bpos - 1).(apos - 1) +. (subst matrix abase bbase) +.
                (go gap aseq apos opening) +. (go gap bseq bpos opening))
            ((!vertical).(bpos - 1).(apos - 1) +. 
                (go' matrix opening gap bbase abase))
            ((!horizontal).(bpos - 1).(apos - 1) +.
                (go' matrix opening gap abase bbase));
    (!horizontal).(bpos).(apos) <-
        min 
            ((!horizontal).(bpos).(apos - 1) +. (ge bbase gap bextension))
            ((!diagonal).(bpos).(apos - 1) +. (ge bbase gap bextension) +.
                (go gap bseq bpos opening));
    (!vertical).(bpos).(apos) <-
        min 
            ((!vertical).(bpos - 1).(apos) +. (ge abase gap aextension))
            ((!diagonal).(bpos - 1).(apos) +. (ge abase gap aextension) +.
                (go gap aseq apos opening));
    (!diagonal).(bpos).(apos) <-
        (diag gap abase bbase) +. 
            (min 
                (!diagonal).(bpos - 1).(apos - 1)
                ((!cost).(bpos - 1).(apos - 1) +. (go gap aseq apos opening) +.
                    (go gap bseq bpos opening)));
    select_best_direction apos bpos

let initialize_aff gap opening matrix aseq bseq =
    let a = 1
    and h = 2
    and v = 4 in
    (!cost).(0).(0) <- 0.;
    (!diagonal).(0).(0) <- large_float;
    (!vertical).(0).(0) <- go gap aseq 1 opening;
    (!horizontal).(0).(0) <- go gap bseq 1 opening;
    (!direction).(0).(0) <- a;
    let alen = Sequence.length aseq in
    for i = 1 to alen - 1 do
        let bbase = Sequence.get bseq i in
        let extension = Cost_matrix.Two_D.cost bbase gap matrix.matrix in
        (!cost).(0).(i) <- large_float;
        (!diagonal).(0).(i) <- large_float;
        (!vertical).(0).(i) <- large_float;
        (!horizontal).(0).(i) <- 
            (!horizontal).(0).(i - 1) +. (ge bbase gap extension);
        (!direction).(0).(i) <- h;
    done;
    let blen = Sequence.length bseq in
    for i = 1 to blen - 1 do
        let abase = Sequence.get aseq i in
        let extension = Cost_matrix.Two_D.cost abase gap matrix.matrix in
        (!cost).(i).(0) <- large_float;
        (!diagonal).(i).(0) <- large_float;
        (!horizontal).(i).(0) <- large_float;
        (!vertical).(i).(0) <- 
            (!vertical).(i - 1).(0) +. (ge abase gap extension);
        (!direction).(i).(0) <- v;
    done;
    ()

let align matrix aseq apos bseq bpos =
    let abase = Sequence.get aseq apos
    and bbase = Sequence.get bseq bpos 
    and gap = Cost_matrix.Two_D.gap matrix.matrix in
    let subs = float_of_int (Cost_matrix.Two_D.cost abase bbase matrix.matrix)
    and adel = float_of_int (Cost_matrix.Two_D.cost abase gap matrix.matrix)
    and bins = float_of_int (Cost_matrix.Two_D.cost gap bbase matrix.matrix) in
    let tot_subs = 
        subs +. (eventf_ac subs (apos - 1) (bpos - 1) matrix) +. 
        (costf matrix (apos - 1) (bpos - 1))
    and tot_adel = 
        adel +. (eventf_ac adel (apos - 1) bpos matrix) +. 
        (costf matrix (apos - 1) bpos)
    and tot_bins =
        bins +. (eventf_ac bins apos (bpos - 1) matrix) +. 
        (costf matrix apos (bpos - 1))
    in
    if tot_subs < tot_adel then begin
        if tot_subs < tot_bins then begin
            (!event).(bpos).(apos) <- 
                matrix.event_cost +.
                (eventf subs (apos - 1) (bpos - 1) matrix);
            (!cost).(bpos).(apos) <- tot_subs;
            (!direction).(bpos).(apos) <- 1;
        end else begin
            (!event).(bpos).(apos) <- 
                matrix.event_cost +.
                (eventf bins apos (bpos - 1) matrix);
            (!cost).(bpos).(apos) <- tot_bins;
            (!direction).(bpos).(apos) <- 2;
        end
    end else if tot_adel < tot_bins then begin
        (!event).(bpos).(apos) <- 
            matrix.event_cost +.
            eventf adel (apos - 1) bpos matrix;
        (!cost).(bpos).(apos) <- tot_adel;
        (!direction).(bpos).(apos) <- 4;
    end else begin 
        (!event).(bpos).(apos) <- 
            matrix.event_cost +.
            eventf bins apos (bpos - 1) matrix;
        (!cost).(bpos).(apos) <- tot_bins;
        (!direction).(bpos).(apos) <- 2;
    end

let align matrix aseq bseq =
    let alen = Sequence.length aseq
    and blen = Sequence.length bseq in
    for i = 1 to alen - 1 do
        for j = 1 to blen - 1 do
            align matrix aseq i bseq j;
        done;
    done

let affine_align matrix aseq bseq =
    let alen = Sequence.length aseq 
    and blen = Sequence.length bseq in
    let opening = 
        float_of_int (Cost_matrix.Two_D.gap_opening matrix.matrix)
    and gap = Cost_matrix.Two_D.gap matrix.matrix in
    for i = 1 to alen - 1 do
        for j = 1 to blen - 1 do
            affine_align gap opening matrix aseq i bseq j;
        done;
    done


let initialize matrix aseq bseq = 
    let alen = Sequence.length aseq
    and blen = Sequence.length bseq 
    and gap = Cost_matrix.Two_D.gap matrix.matrix in
    (!cost).(0).(0) <- 0.;
    (!event).(0).(0) <- matrix.event_cost +. matrix.first_event_cost;
    (!direction).(0).(0) <- 1;
    for i = 1 to blen - 1 do
        let mcost = 
            float_of_int 
            (Cost_matrix.Two_D.cost gap (Sequence.get bseq i) matrix.matrix) 
        in
        (!cost).(i).(0) <- (!cost).(i - 1).(0) +. mcost +.
            (if mcost <> 0. then (!event).(i - 1).(0) else 0.);
        (!direction).(i).(0) <- 2;
        (!event).(i).(0) <- matrix.event_cost +. matrix.first_event_cost;
    done;
    let cost_row = (!cost).(0) 
    and direction_row = (!direction).(0)
    and event_row = (!event).(0) in
    for j = 1 to alen - 1 do
        let cost = 
            float_of_int 
            (Cost_matrix.Two_D.cost gap (Sequence.get aseq j) matrix.matrix) 
        in
        cost_row.(j) <- cost_row.(j - 1) +. cost +. 
            (if cost <> 0. then event_row.(j - 1) else 0.);
        direction_row.(j) <- 4;
        event_row.(j) <- matrix.event_cost +. matrix.first_event_cost;
    done

let allocate_memory msize matrix =
    if msize < !size then ()
    else begin
        cost := Array.make_matrix msize msize 0.;
        vertical := Array.make_matrix msize msize 0.;
        horizontal := Array.make_matrix msize msize 0.;
        diagonal := Array.make_matrix msize msize 0.;
        direction := Array.make_matrix msize msize 0;
        event := Array.make_matrix msize msize 0.; 
        size := msize;
    end

let prepend_gap gap median =
    if gap = Sequence.get median 0 then ()
    else Sequence.prepend median gap

let backtrace matrix aseq bseq =
    let alen = Sequence.length aseq
    and blen = Sequence.length bseq in
    let median = Sequence.create (alen + blen + 2) 
    and a_algn = Sequence.create (alen + blen + 2) 
    and b_algn = Sequence.create (alen + blen + 2) in
    let apos = ref (alen - 1)
    and bpos = ref (blen - 1) in
    let direction = !direction 
    and gap = Cost_matrix.Two_D.gap matrix.matrix in
    while (!apos > 0) && (!bpos > 0) do
        let abase = Sequence.get aseq !apos
        and bbase = Sequence.get bseq !bpos in
        let base_median = Cost_matrix.Two_D.median abase bbase matrix.matrix in
        if base_median <> gap then Sequence.prepend median base_median;
        match direction.(!bpos).(!apos) with
        | 1 -> 
                Sequence.prepend a_algn abase;
                Sequence.prepend b_algn bbase;
                decr apos;
                decr bpos;
        | 2 -> 
                Sequence.prepend a_algn gap;
                Sequence.prepend b_algn bbase;
                decr bpos;
        | 4 ->
                Sequence.prepend a_algn abase;
                Sequence.prepend b_algn gap;
                decr apos;
        | _ -> assert false
    done;
    while !apos > 0 do
        let abase = Sequence.get aseq !apos in
        let base_median = Cost_matrix.Two_D.median gap abase matrix.matrix in
        if base_median <> gap then Sequence.prepend median base_median;
        decr apos;
    done;
    while !bpos > 0 do
        let bbase = Sequence.get bseq !bpos in
        let base_median = Cost_matrix.Two_D.median gap bbase matrix.matrix in
        if base_median <> gap then Sequence.prepend median base_median;
        decr bpos;
    done;
    prepend_gap gap median;
    ((!cost).(blen - 1).(alen - 1)), a_algn, b_algn, median


let align aseq bseq matrix =
    let () = 
        try
            allocate_memory (max (Sequence.length aseq) (Sequence.length bseq)) matrix
        with
        | err -> failwith "1"
    in
    let () = 
        try
    initialize matrix aseq bseq;
        with
        | err -> failwith "2"
    in
    let () = 
        try
    align matrix aseq bseq;
        with
        | err -> failwith "3"
    in
        try
    backtrace matrix aseq bseq 
        with
        | err -> failwith "4"
