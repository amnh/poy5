let cost = ref [||]
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

let align matrix aseq apos bseq bpos =
    let abase = Sequence.get aseq apos
    and bbase = Sequence.get bseq bpos 
    and gap = Cost_matrix.Two_D.gap matrix.matrix in
    let subs = float_of_int (Cost_matrix.Two_D.cost abase bbase matrix.matrix)
    and ains = float_of_int (Cost_matrix.Two_D.cost gap abase matrix.matrix)
    and bins = float_of_int (Cost_matrix.Two_D.cost gap bbase matrix.matrix) in
    let tot_subs = 
        subs +. (eventf_ac subs (apos - 1) (bpos - 1) matrix) +. 
        (costf matrix (apos - 1) (bpos - 1))
    and tot_ains = 
        ains +. (eventf_ac ains (apos - 1) bpos matrix) +. 
        (costf matrix (apos - 1) bpos)
    and tot_bins =
        bins +. (eventf_ac bins apos (bpos - 1) matrix) +. 
        (costf matrix apos (bpos - 1))
    in
    if tot_subs < tot_ains then begin
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
    end else if tot_ains < tot_bins then begin
        (!event).(bpos).(apos) <- 
            matrix.event_cost +.
            eventf ains (apos - 1) bpos matrix;
        (!cost).(bpos).(apos) <- tot_ains;
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
    (!cost).(blen - 1).(alen - 1), a_algn, b_algn, median


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
