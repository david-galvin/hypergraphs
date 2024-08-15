#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(unused_mut)]

use rustc_hash::{FxHashMap, FxHashSet};
use std::fmt;
use std::env;
use std::time::{Duration, Instant};
//use smallvec::{smallvec, SmallVec}; // https://docs.rs/smallvec/1.10.0/smallvec/struct.SmallVec.html


// f_r(t, n) is the minimum number of maximal monochromatic cliques
// across all t-colorings of the complete r-uniform hypergraph on n vertices

// TODO Limit work when changing an edge color (instead of recalculating from scratch)
// TODO Parallelization
// TODO Experiment with annealing -- retain the best hypergraphs, and
//      randomizing edges with slowly decreasing probability
// TODO Grow a graph by:
//        a) cloning a vertex, and only randomizing the other new edges (having the cloned vertex and its clone as two nodes of the edge)
//        b) creating a new vertex and only randomizing new edge colors.
// TODO Develop a likely-unique signature for each graph, based on its degree & clique sequence be color (count of degrees, count of cliques)
// println!("Time elapsed: {:?}", start.elapsed());

fn main() {
  let start = Instant::now();
  // Get user inputs. Since we use f_{edge order}(color_count, graph_order) in the paper, we'll preserve this ordering
  let args: Vec<String> = env::args().collect();
  if args.len() < 3 {
    println!("\ncargo run edge_order color_ct graph_order. E.g., f_3(2, 6) = cargo run 3 2 6\n");
    return;
  }
  
  // Inputs:
  // Edge Order: The number of vertices that define an edge. 2 for a normal graph.
  // Color Count: Each edge has one of color_ct colors
  // Graph Order: The number of vertices in the graph
  let edge_order: u8 = args[1].parse().unwrap();
  let color_ct: u8 = args[2].parse().unwrap();
  let graph_order: u8 = args[3].parse().unwrap();
  
  let randomize_colors: bool;
  
  // Derived:
  // Graph Size: The count of edges in the graph
  let graph_size: usize = math::choose(graph_order as usize, edge_order as usize);

  // Calculate & print upper bound on count of maximal monochromatic cliques
  let (upper_bound, theorem) = get_upper_bound(graph_order, edge_order, color_ct);
  
  let threshold_str = format!("{}", theorem);
  println!("{}", threshold_str);
  
  let mut h: HyperGraph = HyperGraph::new(edge_order, color_ct, graph_order);
  
  
  // set edge colors if they've been supplied
  // idea: we'll color edges from left to right, and move any intentionally 
  // uncolored ones to the right. If at the end there are any uncolored ones
  // left, those will be randomly explored.
  let mut arg_index: usize = 4;
  let mut edge_index_right: usize = h.graph_size; // refers to the left-most edge with color
  let mut edge_index_left: usize = h.graph_size; // refers to the left-most edge not yet checked for color
  let mut color: u8;
  while args.len() > arg_index {
    color = (arg_index - 4) as u8;
    edge_index_left = 0;
    'edge_loop: while edge_index_left < edge_index_right {
      let edge = &mut h.edges[edge_index_left];
      for clique_members in args[arg_index].parse::<String>().unwrap().split(",") {
        let clique_members_u64: u64 = clique_members.parse().expect("Not a valid number");
        if edge.members & clique_members_u64 == edge.members {
          edge.color = color;
          edge_index_right -= 1;
          h.edges.swap(edge_index_left, edge_index_right);
          continue 'edge_loop;
        }
      }
      edge_index_left += 1;
    }
    arg_index += 1;
  }
  
  if edge_index_left == 0 {
    randomize_colors = false;
    h.find_cliques_from_scratch();
    println!("\n\nInput Hypergraph had all edge colors specified:");
    println!("{}", h);
  } else {
    h.last_random_edge_index = edge_index_left - 1;
    randomize_colors = true;
  }
  
  // Show all hypergraphs with uniform coloring of all unspecified edges
  if randomize_colors && h.last_random_edge_index < graph_size - 1 {
    for cur_color in 0..color_ct {
      for i in 0..=h.last_random_edge_index {
        h.edges[i].color = cur_color;
      }
      h.find_cliques_from_scratch();
      println!("\n\nInitial Hypergraph with unspecified edges colored {}:", cur_color);
      println!("{}", h);
    }
  }
  
  
  let mut loops_without_improvement: usize;
  let mut edge_index_to_try: usize = 0;
  let mut best_from_current_start: usize;
  let mut best_from_all_starts: usize = usize::MAX;
  let mut annealing_count: usize = 0;
  
  loop {
    if randomize_colors {
      h.randomize_colors();
    } else {
      break;
    }
    h.find_cliques_from_scratch();
    best_from_current_start = h.maximal_color_clique_ct;

    loops_without_improvement = 0;
    
    while loops_without_improvement < graph_size  {
      //println!("loops without improvement: {}", loops_without_improvement);
      loops_without_improvement += 1;

      // Try all alternate colors for the current edge
      for _ in 0..(color_ct - 1) {
        h.increment_edge_color(edge_index_to_try);
        h.find_cliques_from_scratch();
        if h.maximal_color_clique_ct < best_from_current_start {
          if h.maximal_color_clique_ct < best_from_all_starts || 
          ((h.maximal_color_clique_ct == best_from_all_starts) && h.maximal_color_clique_ct <= upper_bound) {
            
            println!(
              "\n------------------------------------------\n\n{}.\n\nIMPROVEMENT: {} -> {} (flipped edge {}). Annealings: {}, Time: {:?}", 
              threshold_str, 
              best_from_all_starts, 
              h.maximal_color_clique_ct,
              edge_index_to_try, 
              annealing_count,
              start.elapsed());
            
            best_from_all_starts = h.maximal_color_clique_ct;
            println!("{}", h);
          }
          loops_without_improvement = 0;
          best_from_current_start = h.maximal_color_clique_ct
        }
      }
      
      if loops_without_improvement > 0 {
        // No improvement found; revert the color
        h.increment_edge_color(edge_index_to_try);
      }
      
      // Shift to the next edge index
      edge_index_to_try = (edge_index_to_try + 1) % graph_size;
    }
    annealing_count += 1;
    
    if (annealing_count <= 100) || 
       (annealing_count <= 1000 && annealing_count % 10 == 0) ||
       (annealing_count <= 10000 && annealing_count % 100 == 0) ||
       (annealing_count <= 100000 && annealing_count % 1000 == 0) ||
       (annealing_count <= 1000000 && annealing_count % 10000 == 0) ||
       (annealing_count <= 10000000 && annealing_count % 100000 == 0) ||
       (annealing_count <= 100000000 && annealing_count % 1000000 == 0) ||
       (annealing_count >  100000000 && annealing_count % 10000000 == 0) {
      println!("annealings: {}, Time: {:?}", annealing_count, start.elapsed());
    }
  }
}

struct HyperGraph {
  edges: Vec<Clique>,
  cliques: FxHashMap<u8, Vec<Clique>>,
  members_of_cliques_which_should_be_deactivated: FxHashSet<u64>,
  maximal_color_clique_ct: usize,
  edge_order: u8,
  color_ct: u8,
  graph_order: u8,
  graph_size: usize,
  last_random_edge_index: usize,
  set_bit_getter: math::SetBitGetter,
}

impl HyperGraph {
  fn new(edge_order: u8, color_ct: u8, graph_order: u8) -> HyperGraph {

    let cliques = FxHashMap::<u8, Vec<Clique>>::default();
    
    let graph_size: usize = math::choose(graph_order as usize, edge_order as usize);
  
    let h = HyperGraph {
      edges: Clique::generate_all_cliques(u8::max_value(), edge_order, graph_order),
      cliques,
      members_of_cliques_which_should_be_deactivated: FxHashSet::<u64>::default(),
      maximal_color_clique_ct: 0,
      edge_order,
      color_ct,
      graph_order,
      graph_size: graph_size,
      last_random_edge_index: graph_size - 1, 
      set_bit_getter: math::SetBitGetter::new(),
    };

    h
  }
  
  fn get_key(&self, color: u8, order: u8) -> u8 {
    color * self.graph_order + order
  }

  fn randomize_colors(&mut self) {
    for i in 0..=self.last_random_edge_index {
      self.edges[i].set_color(fastrand::u8(0..self.color_ct));
    }
  }
  
  fn random_clique_growth(&mut self) {
	// TODO Ensure clique doesn't contain all vertices.
	
	// TODO Pick random vertex  vnot in clique
	
	// TODO color all edges between every set of r-1 of the clique's vertices & v
	//      with the clique's color.
	  
	  
	// cliques: FxHashMap<u8, Vec<Clique>>,
	let mut cliques_count: usize = 0;
	let mut growth_clique: &Clique;
	let target: usize = ((fastrand::f64()) * (self.maximal_color_clique_ct as f64)) as usize;
	for (key, val) in self.cliques.iter() {
      cliques_count += val.len();
      if cliques_count >= target {
		growth_clique = &val[cliques_count % target];
		break;
	  }
    }
  }
  
  fn increment_edge_color(&mut self, edge_index: usize) {
    self.edges[edge_index].increment_color(self.color_ct);
  }
  
  fn confirm_cliques_vec_initialized(&mut self, key: u8) {
    self.cliques.entry(key).or_insert_with(Vec::<Clique>::new);
  }
  
  fn add_clique(&mut self, clique: Clique) {
    self.maximal_color_clique_ct += 1;
    let key = self.get_key(clique.color, clique.order);
    self.confirm_cliques_vec_initialized(key);
    self.cliques.get_mut(&key).expect("Uninitalized Vector").push(clique);
  }
  
  fn add_clique_vec(&mut self, mut clique_vec: Vec<Clique>, key: u8) {
    self.maximal_color_clique_ct += clique_vec.len();
    self.confirm_cliques_vec_initialized(key);
    self.cliques.get_mut(&key).expect("Uninitalized Vector").append(&mut clique_vec);
  }
  
  fn mark_nonmaximal_cliques_inactive(&mut self, color: u8, order_of_larger_clique: u8) {
    let key_big = self.get_key(color, order_of_larger_clique);
    if !&self.cliques.contains_key(&key_big) {
      return;
    }
    
    let key_small = self.get_key(color, order_of_larger_clique - 1);
    
    for big_clique in &self.cliques[&key_big] {
      self.set_bit_getter.reset(big_clique.members);
      while self.set_bit_getter.are_there_more_bits() {
        self.members_of_cliques_which_should_be_deactivated.insert(big_clique.members ^ self.set_bit_getter.get_bit());
      }
    }


    // Parse all small cliques, deleting any whose members are marked for deletion.
    for small_clique in self.cliques.get_mut(&key_small).expect("Uninitalized Vector") {
      if small_clique.is_active && self.members_of_cliques_which_should_be_deactivated.contains(&small_clique.members) {
        self.members_of_cliques_which_should_be_deactivated.remove(&small_clique.members);
        small_clique.set_inactive();
        self.maximal_color_clique_ct -= 1;
      }
    }
  }


/*
  fn find_cliques_after_edge_recoloring(&mut self, pri_color: u8, edge_index: usize) {
    let mut key: u8;
    let mut key_small: u8;
    let mut any_activated: bool;
    let edge_members = self.edges[edge_index].members;
    let mut buffer: Vec<_> = Vec::with_capacity(128); // Preallocate buffer with a reasonable size
    
    for order in self.edge_order..=self.graph_order {
      key = self.get_key(pri_color, order);
      key_small = self.get_key(pri_color, order - 1);
      any_activated = false;

      if let Some(cliques) = self.cliques.get_mut(&key) {
        buffer.clear(); // Reuse buffer
        
        // CHECK FOR ACTIVE CLIQUES CONTAINING THE EDGE AS IF IT WERE STILL ITS OLD/PRIOR COLOR
        // STORE IN A BUFFER FOR DEACTIVIATION
        for clique in cliques.iter_mut() {
          if clique.is_active && (edge_members & clique.members) == edge_members {
            clique.set_inactive();
            buffer.push(clique.members); // Store members for the second loop
          }
        }
        
        // ACTIVATE ANCESTOR CLIQUES
        if let Some(small_cliques) = self.cliques.get_mut(&key_small) {
          for pri_color_small_clique in small_cliques.iter_mut() {
            if buffer.iter().any(|&members| (pri_color_small_clique.members & members) == pri_color_small_clique.members) {
              pri_color_small_clique.set_active();
              any_activated = true;
            }
          }
        }
      }

      // CONFIRM ACTIVATED ANCESTOR CLIQUES AREN'T CONTAINED IN LARGER ACTIVE CLIQUES
      if any_activated {
        self.mark_nonmaximal_cliques_inactive(pri_color, order);
      }
    }
  }

*/
  
  fn find_cliques_from_scratch(&mut self) {
    self.cliques.clear();
    self.maximal_color_clique_ct = 0;
    let mut smaller_cliques_count: usize;
    let mut order_usize: usize;
    let mut compatible_cliques_needed: usize;
    let mut key_small_cliques: u8;  // key for order-1 cliques
    let mut explored_subsets = FxHashSet::<u64>::default();
    
    let mut cliques_to_add = Vec::<Clique>::new();

    for color in 0..self.color_ct {
      // println!("\n  color: {}", color);
      
      // Create all trivial color-cliques (those having order equal to EDGE_ORDER - 1)
      self.add_clique_vec(
        Clique::generate_all_cliques(color, self.edge_order - 1, self.graph_order), 
        self.get_key(color, self.edge_order - 1));

      // Create all color-edge-cliques
      let mut edge_cliques = Vec::<Clique>::new();
      for edge in &mut self.edges {
        if edge.color == color {
          edge_cliques.push(edge.clone());
        }
      }
      self.add_clique_vec(
        edge_cliques,
        self.get_key(color, self.edge_order));
        
      
      // Mark fully contained trivial color-cliques as inactive (aka non-maximal)  
      self.mark_nonmaximal_cliques_inactive(color, self.edge_order);
      
      // Create all color-order-cliques for orders > edge_order
       
      // 'smaller_cliques' refers to all cliques that match the color
      // we're examining, and have (order-1) vertices.
      // The k'th bit of a cliques_index corresponds to the index in cliques[(color, order - 1)]
      
      // These are used in a hot loop; allocating memory once here.
      let mut candidate_clique_members: u64;
      let mut cur_smaller_clique_members: u64;
      let mut i_members: u64;
      let mut j_members: u64;
      
      for order in (self.edge_order + 1)..(self.graph_order + 1) {
        // println!("    order: {}", order); 
        
        order_usize = order as usize;
        key_small_cliques = self.get_key(color, order - 1);
        
        // Confirm we have enough smaller cliques to build a bigger clique out of
        let smaller_cliques = match self.cliques.get(&key_small_cliques) {
          Some(cliques) => cliques,
          None => break,
        };
        
        smaller_cliques_count = smaller_cliques.len();
        if smaller_cliques_count < order_usize {
          break;
        }

        cliques_to_add.clear(); // Clear previous cliques
        explored_subsets.clear(); // Clear explored subsets

        for i in 0..=(smaller_cliques_count - order_usize) {
          i_members = smaller_cliques[i].members;
          for j in (i + 1)..=(smaller_cliques_count - order_usize + 1) {
            j_members = smaller_cliques[j].members;
            candidate_clique_members = i_members | j_members;
            if !(j_members ^ candidate_clique_members).is_power_of_two() || explored_subsets.contains(&candidate_clique_members) {
              continue;
            }
            explored_subsets.insert(candidate_clique_members);
            compatible_cliques_needed = order_usize - 2;
            for k in (j + 1)..smaller_cliques_count {
              cur_smaller_clique_members = smaller_cliques[k].members;
              if cur_smaller_clique_members & candidate_clique_members == cur_smaller_clique_members {
                compatible_cliques_needed -= 1;
                if compatible_cliques_needed == 0 {
                  // We've found an order-clique!
                  cliques_to_add.push(Clique::new(candidate_clique_members, color, self.graph_order));
                }
                continue;
              }
            }
          }
        }

        // Move cliques_to_add into add_clique
        while let Some(clique) = cliques_to_add.pop() {
          self.add_clique(clique);
        }


        // Having finished finding color-order cliques, mark any nonmaximal smaller cliques as inactive
        self.mark_nonmaximal_cliques_inactive(color, order);
      }
    }
  }

/*
  fn reset_cliques_after_edge_recoloring(&mut self) {
    for color in 0..self.color_ct {
      for clique in self.cliques.get_mut(&self.get_key(color, self.edge_order - 1)).expect("Missing clique vector") {
        clique.set_active();
      }
      
      for order in (self.edge_order + 1)..(self.graph_order + 1) {
        self.cliques.get_mut(&self.get_key(color, order)).expect("Missing clique vector").clear();
      }
    }
  }
*/

  fn get_string(&self) -> String {
    let mut vertex_clique_counts: Vec<u8> = vec![0; (self.color_ct * self.graph_order) as usize];
    let mut cliques_str = format!("cargo run --release {} {} {}", self.edge_order, self.color_ct, self.graph_order);
    let mut ret_str = format!(
      "\nH: A complete {}-uniform hypergraph on {} vertices with {} colors, cliques: {}, edges: {}\n", 
      self.edge_order, 
      self.graph_order, 
      self.color_ct,
      self.maximal_color_clique_ct,
      self.graph_size);
      
    // Printing edges is too noisy!  
    /*for edge in &self.edges {
      ret_str += &format!("{}\n", edge);
    }*/

    ret_str += "\nMaximal Color-Cliques\n";
    for color in 0..self.color_ct {
      ret_str += &format!("  color {}:\n", color);
      cliques_str += &format!("");
      let mut first_of_color: bool = true;
      for order in (self.edge_order - 1)..(self.graph_order + 1) {
        let mut order_str = format!("    order {}:\n", order);
        let mut include_order_str: bool = false;
        if !&self.cliques.contains_key(&self.get_key(color, order)) {
          continue;
        }
        for clique in &self.cliques[&self.get_key(color, order)] {
          if !clique.is_active {
            continue;
          }
          if order >= self.edge_order {
            if first_of_color {
              cliques_str += &format!(" ");
              first_of_color = false;
            } else {
              cliques_str += &format!(",");
            }
            cliques_str += &format!("{}", clique.members);
            
          }
          
          include_order_str = true;
          order_str += &format!("      {}\n", clique);
          for i in 0..self.graph_order {
            if (1 << i) & clique.members != 0 {
              vertex_clique_counts[(color * self.graph_order + i) as usize] += 1;
            }
          }
        }
        if include_order_str {
          ret_str += &order_str;
        }
      }
    }
    ret_str += &format!("\n  {} Maximal Color Cliques Found. Vertices on few cliques: (", self.maximal_color_clique_ct);
    for color in 0..self.color_ct {
      let mut min_count = self.graph_order;
      for i in (color * self.graph_order)..((color + 1) * self.graph_order) {
        if vertex_clique_counts[i as usize] < min_count {
          min_count = vertex_clique_counts[i as usize];
        }
      }
      ret_str += &format!("{}", min_count);
      if color != self.color_ct - 1 {
        ret_str += ", ";
      }
    }
    ret_str += ")\n";
    ret_str += &cliques_str;
    ret_str
  }
}

impl fmt::Display for HyperGraph {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{}", self.get_string())
  }
}

struct Clique {
  members: u64,
  color: u8,
  order: u8,
  graph_order: u8,
  is_active: bool,
}

impl Clique {
  fn new(members: u64, color: u8, graph_order: u8) -> Clique {
    Clique {
      members,
      color,
      order: members.count_ones() as u8,
      graph_order,
      is_active: true,
    }
  }
  
  fn generate_all_cliques(clique_color: u8, clique_order: u8, graph_order: u8) -> Vec<Clique> {
    let mut ret_vec: Vec<Clique> = Vec::new();
    let mut cur_members: u64 = (1 << clique_order) - 1;
    let max_members: u64 = (1 << graph_order) - (1 << (graph_order - clique_order));
    while cur_members <= max_members {
      ret_vec.push(Clique::new(cur_members, clique_color, graph_order));
      math::set_next_bit_permutation_u64(&mut cur_members);
    }
    ret_vec
  }
  
  fn set_color(&mut self, new_color: u8) {
    self.color = new_color;
  }
  
  fn increment_color(&mut self, color_ct: u8) {
    self.color = (self.color + 1) % color_ct;
  }
  
  fn set_active(&mut self) {
    self.is_active = true;
  }
  
  fn set_inactive(&mut self) {
    self.is_active = false;
  }
  
  fn get_string(&self) -> String {
    let mut ret_str = String::new();
    for i in (0..self.graph_order).rev() {
      if self.members & (1 << i) > 0 {
        ret_str += "\u{25AA}";
      } else {
        ret_str += "\u{2B1D}";
      }
    }
    ret_str += " ";
    let members_str = format!("{0:b}", self.members);
    ret_str += format!("{:0>width$}", members_str, width = self.graph_order as usize).as_str();
    ret_str += format!(" ({})", self.members).as_str();
    ret_str += format!(" ({})", self.color).as_str();
    if !self.is_active {
      ret_str += " (I)"; // inactive
    }
    ret_str
  }
  
  fn clone(&self) -> Clique {
    Clique::new(self.members, self.color, self.graph_order)
  }
}

impl fmt::Display for Clique {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.get_string())
    }
}


// f_r(t, n)
// get the upper bound of the minimum amount of maximal monochromatic cliques based on the Galvin / Galvin / Erdos / Krieger paper
fn get_upper_bound(n: u8, r: u8, t: u8) -> (usize, String) { // (bound, theorem, is_equal (vs upper bound))
  let n = n as usize;
  let r = r as usize;
  let t = t as usize;
  
  let mut upper_bound: usize = usize::MAX;
  let mut cur_upper_bound: usize;
  let mut theorem_str = String::new();
  
  // Theorem 1
  if r == 2 && t == 2 {
    cur_upper_bound = n + 1;
    if cur_upper_bound < upper_bound {
      upper_bound = cur_upper_bound; 
      theorem_str = format!("Thm 1: f_2(2, {}) = {}\n       f_2(2, n) = n + 1", n, cur_upper_bound);
    }
  }
  
  // Theorem 2
  if r >= 2 && t >= 2 {
    
    // 2.a
    if 1 <= n && n < r {
      cur_upper_bound = t;
      if cur_upper_bound < upper_bound {
        upper_bound = cur_upper_bound; 
        theorem_str = format!("Thm 2a: f_{}({}, {}) = {}\n        f_r(t, n) = t if 1 <= n < r", r, t, n, cur_upper_bound);
      }
    }
    
    // 2.b
    if n == r {
      cur_upper_bound = (t-1) * n + 1;
      if cur_upper_bound < upper_bound {
        upper_bound = cur_upper_bound; 
        theorem_str = format!("Thm 2b: f_{}({}, {}) = {}\n        f_r(t, n) = (t-1) * n + 1 if n == r", r, t, n, cur_upper_bound);
      }
    }
    
    // 2.c NOT IMPLEMENTED
    if n == r + 1 {
      let turan_edge_count: usize = ((n / t).pow((t - n.rem_euclid(t)) as u32)) * ((n / t + 1).pow(n.rem_euclid(t) as u32));
      cur_upper_bound = math::choose(n + 1, 2) + (t - 2) * math::choose(n, 2) - turan_edge_count;
      if cur_upper_bound < upper_bound {
        upper_bound = cur_upper_bound; 
        theorem_str = format!("Thm 2c: f_{}({}, {}) = {}\n        f_r(t, n) = choose(n+1, 2) + (t-2) * choose(n, 2) - e(T_n,t)", r, t, n, cur_upper_bound);
      }
    }
    
    // 2.d
    if n == r + 1 && t == 2 {
      cur_upper_bound = ((n + 1).pow(2)) / 4;
      if cur_upper_bound < upper_bound {
        upper_bound = cur_upper_bound; 
        theorem_str = format!("Thm 2d: f_{}({}, {}) = {}\n        f_r(t, n) = floor((n + 1)^2 / 4) if n == r + 1 and t == 2", r, t, n, cur_upper_bound);
      }
    }
  }
  
  // Theorem 6
  if r >= 2 && t >= 2 && n >= r-1 {
    cur_upper_bound = (t - 1) * math::choose(n, r - 1) + 1;
    if cur_upper_bound < upper_bound {
     upper_bound = cur_upper_bound; 
     theorem_str = format!("Thm 6: f_{}({}, {}) <= {}\n       f_r(t, n) <= (t-1) * choose(n, r-1) + 1 if r, t >= 2 and n >= r-1", r, t, n, cur_upper_bound);
    }
  }
  
  
  if n >= 2 && r == 3 && t == 2 {
    // Theorem 8
    cur_upper_bound = (n + 1).pow(2) / 4;
    if cur_upper_bound < upper_bound {
     upper_bound = cur_upper_bound; 
     theorem_str = format!("Thm 8: f_{}({}, {}) <= {}\n       f_3(2, n) <= floor((n + 1)^2 / 4)) if n >= 2", r, t, n, cur_upper_bound);
    }
    
    // Lemma 9
    if n >= 7 {
      cur_upper_bound = n.pow(2)/4 + 2;
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Lem 9: f_{}({}, {}) <= {}\n       f_3(2, n) <= floor(n^2 / 4) + 2 if n >= 7", r, t, n, cur_upper_bound);
      }
    }
    
    // Theorem 11
    if n >= 7 {
      cur_upper_bound = (n + 1).pow(2) / 4 - 2;
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Thm 11: f_{}({}, {}) <= {}\n        f_3(2, n) <= floor((n + 1)^2 / 4) - 2 if n >= 7", r, t, n, cur_upper_bound);
      }
    }
    
    // Theorem 12
    if n >= 15 {
      cur_upper_bound = n.pow(2) / 4 + 5;
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Thm 12: f_{}({}, {}) <= {}\n        f_3(2, n) <= floor(n^2 / 4) + 5", r, t, n, cur_upper_bound);
      }
    }
  }
    
  if r == 2 {
    // Theorem 13
    if t >= 2 && n >= 1 {
      cur_upper_bound = (t - 1) * n + 1;
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Thm 13: f_{}({}, {}) <= {}\n        f_2(t, n) <= (t-1) * n + 1 if t >= 2 and n >= 1", r, t, n, cur_upper_bound);
      }
    }
    
    // Theorem 14
    if n <= 2 * (t/2) + (t % 2) {
      cur_upper_bound = t * n - math::choose(n, 2);
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Thm 14: f_{}({}, {}) = {}\n        f_2(t, n) = t * n - choose(n, 2) if n <= 2 * ceiling(t/2)", r, t, n, cur_upper_bound);
      }
    }
    
    // Corollary 15a
    if t == n + 1 {
      cur_upper_bound = n * (n + 3) / 2;
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Cor 15a: f_{}({}, {}) = {}\n         f_2(n + 1, n) = n(n + 3)/2", r, t, n, cur_upper_bound);
      }  
    }
    
    // Corollary 15b
    if (t % 2 == 1) && ((t == n-1) || (t == n)) {
      cur_upper_bound = math::choose(t + 1, 2);
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Cor 15b: f_{}({}, {}) = {}\n         f_2(n, n + 1) = choose(n + 1, 2) = f_2(n, n) if n is odd", r, t, n, cur_upper_bound);
      }      
    }
    
    // Corollary 15c
    if (t % 2 == 0) && (t == n) {
      cur_upper_bound = math::choose(n + 1, 2);
      if cur_upper_bound < upper_bound {
       upper_bound = cur_upper_bound; 
       theorem_str = format!("Cor 15c: f_{}({}, {}) = {}\n         f_2(n, n) = choose(n + 1, 2) if n is even", r, t, n, cur_upper_bound);
      } 
    }
  }
  
  // ppo = projective plane order
  let ppos = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128,131,137,139,149,151,157,163,167,169,173,179,181,191,193,197,199,211];
  let mut t_is_ppo_plus_one: bool = false;
  for q in ppos {
    if t == q + 1 {
      t_is_ppo_plus_one = true;
      break;
    }
    if t < q + 1 {
      break;
    }
  }
  
  if r == 2 {
    if t_is_ppo_plus_one {
      let q = t - 1;
      // Theorem 19a
      if (q-1).pow(2) < n && n <= (q-1) * q {
        cur_upper_bound = q.pow(2) + q - 1;
        if cur_upper_bound < upper_bound {
          upper_bound = cur_upper_bound; 
          theorem_str = format!("Thm 19a: f_{}({}, {}) = {}\n         f_2(q + 1, n) = q^2 + q - 1 if there is a PP of order q and (q - 1)^2 < n <= (q - 1)q", r, t, n, cur_upper_bound);
        } 
      }
      // Theorem 19b
      if (q-1) * q < n && n <= q.pow(2) {
        cur_upper_bound = q.pow(2) + q;
        if cur_upper_bound < upper_bound {
          upper_bound = cur_upper_bound; 
          theorem_str = format!("Thm 19b: f_{}({}, {}) = {}\n         f_2(q + 1, n) = q^2 + q if there is a PP of order q and (q-1)q < n <= q^2", r, t, n, cur_upper_bound);
        } 
      }
      // Theorem 20
      if n == (r - 1).pow(2) {
        cur_upper_bound = (r-1).pow(2) + r - 1;
        if cur_upper_bound < upper_bound {
          upper_bound = cur_upper_bound; 
          theorem_str = format!("Thm 20: f_{}({}, {}) = {}\n        f_2(n + 1, n^2) = n^2 + n if there is a PP of order n", r, t, n, cur_upper_bound);
        } 
      }
    }
  }
  
  // Theorem 23
  if r == 2 && t == 3 {
    if n % 3 == 1 {
      cur_upper_bound = n + 2;
    } else {
      cur_upper_bound = n + 3;
    }
    if cur_upper_bound < upper_bound {
      upper_bound = cur_upper_bound; 
      theorem_str = format!("Thm 23: f_{}({}, {}) <= {}\n        f_2(3, n) <= {{ n + 2 if n % 3 == 1\n                     {{ n + 3 otherwise\n        (equality shown for n <= 10)", r, t, n, cur_upper_bound);
    } 
  }

  // Theorem 25
  if r == 2 && t == 4 {
    let n_mod_3 = n % 3;
    if n_mod_3 == 1 {
      cur_upper_bound = n + 3;
    } else if n_mod_3 == 0 {
      cur_upper_bound = n + 4;
    } else if n_mod_3 == 2 || n_mod_3 == 6 || n_mod_3 == 7 {
      cur_upper_bound = n + 5;
    } else {
      cur_upper_bound = n + 6;
    }
    if cur_upper_bound < upper_bound {
      upper_bound = cur_upper_bound; 
      theorem_str = format!("Thm 25: f_{}({}, {}) <= {}\n        f_2(4, n) <= {{ n + 3 if n % 8 == 1\n                     {{ n + 4 if n % 8 == 0\n                     {{ n + 5 if n % 8 == 2,6,7\n                     {{ n + 6 if n % 8 == 3,4,5\n        (equality shown for n <= 10)", r, t, n, cur_upper_bound);
    } 
  }

  // Theorem 27
  if r == 2 && t == 5 && n >= 37{
    let n_mod_5 = n % 5;
    if n_mod_5 == 1 {
      cur_upper_bound = n + 4;
    } else if n_mod_5 == 2 {
      cur_upper_bound = n + 7;
    } else if n_mod_5 == 3 {
      cur_upper_bound = n + 7;
    } else if n_mod_5 == 4 {
      cur_upper_bound = n + 6;
    } else {
      cur_upper_bound = n + 5;
    }
    if cur_upper_bound < upper_bound {
      upper_bound = cur_upper_bound; 
      theorem_str = format!("Thm 27: f_{}({}, {}) = {}\n        f_2(5, n) <= {{ n + 4 if n % 5 == 1\n                     {{ n + 7 if n % 5 == 2,3\n                     {{ n + 6 if n % 5 == 4\n                     {{ n + 5 if n % 5 == 0", r, t, n, cur_upper_bound);
    } 
  }
  (upper_bound, theorem_str)
}
