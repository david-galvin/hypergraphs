#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(unreachable_code)]

use rustc_hash::{FxHashMap, FxHashSet};
use std::fmt;
use std::env;
use std::time::{Instant};
//use smallvec::{smallvec, SmallVec}; // https://docs.rs/smallvec/1.10.0/smallvec/struct.SmallVec.html


// f_r(t, n) is the minimum number of maximal monochromatic cliques
// across all t-colorings of the complete r-uniform hypergraph on n vertices

// TODO Instead of flipping all edges, try having each clique absorbing each vertex not in it
//      and reverting if it doesn't improve.
// TODO Limit work when changing an edge color (instead of recalculating from scratch)
// TODO Parallelization
// TODO Experiment with annealing -- retain the best hypergraphs, and
//      randomizing edges with slowly decreasing probability
// TODO Grow a graph by:
//        a) cloning a vertex, and only randomizing the other new edges (having the cloned vertex and its clone as two nodes of the edge)
//        b) creating a new vertex and only randomizing new edge colors.
//        c) creating a new vertex, find the color with the smallest maximum order, and extend all such color-order cliques to include the new vertex, then randomly color other edges.
// TODO Develop a likely-unique signature for each graph, based on its degree & clique sequence be color (count of degrees, count of cliques)
// TODO Optimize based on lemma9_n1000

fn main() {
  // ARGS 1: GET ARGS
  let args: Vec<String> = env::args().collect();
  if args.len() < 3 {
    println!("\ncargo run edge_order color_count graph_order. E.g., f_3(2, 6) = cargo run 3 2 6\n");
    return;
  }
  
  // ARGS 2: PARSE BASIC ARGS
  let edge_order: u8 = args[1].parse().unwrap();  // The number of vertices that define an edge. Normal graphs have edge_order = 2.
  let color_count: u8 = args[2].parse().unwrap();    // The number of colors, which are integers: 0, 1, ..., color_count - 1.
  let graph_order: u8 = args[3].parse().unwrap(); // The number of vertices in the graph.
  
  // CREATE HYPERGRAPH
  let mut h: HyperGraph = HyperGraph::new(edge_order, color_count, graph_order);
  
  // ARGS 3: PARSE CLIQUE CSV ARGS
  let mut arg_index: usize = 4;
  let mut edge_index_right: usize = h.graph_size; // refers to the left-most edge with color
  let mut edge_index_left: usize = h.graph_size; // refers to the left-most edge not yet checked for color
  let mut color: u8;
  while args.len() > arg_index {
    color = (arg_index - 4) as u8;
    edge_index_left = 0;
    'edge_loop: while edge_index_left < edge_index_right {
      let edge = &mut h.edges[edge_index_left];
      for clique_members in args[arg_index].parse::<String>().unwrap().split(',') {
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
  
  // WERE ANY EDGE COLORS UNSPECIFIED?
  let has_unspecified_edges: bool = edge_index_left > 0;
  
  // IF ALL EDGE COLORS WERE SPECIFIED, PRINT THE GRAPH AND EXIT
  if !has_unspecified_edges {
    h.find_cliques_from_scratch();
    println!("\n\nInput Hypergraph had all edge colors specified:");
    println!("{}", h);
    return;
  } else {
    // OTHERWISE SAVE THE INDEX OF THE LAST UNSPECIFIED EDGE
    h.last_random_edge_index = edge_index_left - 1;
    h.randomize_edge_colors();
    h.find_cliques_from_scratch();
  }


  let mut loops_without_improvement: usize;
  let mut edge_index_to_try: usize = 0;
  let mut key: u8;
  let mut nonmembers: u64;
  let mut vertex_to_add: u64;
  let mut last_key_with_improvement: u8 = 0;
  let mut color: u8;
  let mut order: u8;
  let mut index: usize;
  let mut has_improved: bool;
  let mut new_set_bit_getter = math::SetBitGetter::new();
  let mut skip_reset: bool = false;
  
  // INFINITE LOOP TO LOOK FOR BETTER COLORINGS
  loop {
    h.randomize_edge_colors();
    h.find_cliques_from_scratch();
    h.best_anneal_clique_count = h.current_clique_count;
    if h.current_clique_count < h.best_global_clique_count {
      h.best_global_clique_count = h.current_clique_count;
    }

    key = h.get_key(0, h.edge_order - 1);
    
    loop {
      
      if h.cliques.contains_key(&key) {
        index = 0;
        while index < h.cliques[&key].len() {
          nonmembers = h.cliques[&key][index].members ^ ((1 << h.graph_order) - 1);
          if nonmembers == 0 {
            break;
          }
          
          new_set_bit_getter.reset(nonmembers);
          
          while new_set_bit_getter.are_there_more_bits() {
            has_improved = false;
            vertex_to_add = new_set_bit_getter.get_bit();
            h.grow_clique(key, index, vertex_to_add); 
            h.find_cliques_from_scratch();
            // check for a good result
            if h.current_clique_count <= h.best_anneal_clique_count {
              
              //  WE HAVE IMPROVED GLOBALLY
              if h.current_clique_count < h.best_global_clique_count || 
              ((h.current_clique_count == h.best_global_clique_count) && h.current_clique_count <= h.upper_bound) {
                
                // print a good result
                h.print_header();
                h.best_global_clique_count = h.current_clique_count;
                println!("{}", h);
              }
            }  
                       
            // WE HAVE IMPROVED ON THE BEST RESULT IN THIS ANNEALING
            if h.current_clique_count < h.best_anneal_clique_count {
              last_key_with_improvement = key;
              h.best_anneal_clique_count = h.current_clique_count;
              index = 0;
              has_improved = true;
              if index >= h.cliques[&key].len() {
                break;
              }
            }

            
            // revert if we didn't improve on the current annealing's best
            // TODO: SHOULD THIS BE RANDOM ON EQUALITY?
            // WE HAVE NOT IMPROVED ON THE BEST RESULT IN THIS ANNEALING
            if !has_improved {
              h.revert_clique_growth();
              // TODO FIND A WAY NOT TO find_clique() here!!!
              h.find_cliques_from_scratch();
            }
          }
          index += 1;
        }
      }
      
      // Update color & order; if we've cycled through everything without
      // any improvement, break out so we can start over with a fresh random coloring
      h.rotate_growth_key(&mut key);
      if key == last_key_with_improvement {
        if skip_reset {
          skip_reset = false;
        } else { 
          break;
          //return;
        }
      }
    }

    
    // FOR EVERY CLIQUE
    //   FOR EVERY VERTEX NOT IN THAT CLIQUE
    //     HOW DOES ABSORBING THE VERTEX CHANGE THINGS?
    //       BETTER: RECORD THE (CLIQUE, VERTEX),  AND CONTINUE UNTIL WE REACH IT AGAIN
    //       WORSE: REVERT THE ABSORPTION
    //       UNCH: TBD; MAYBE TREAT AS WORSE, MAYBE RANDOMLY FOLLOW ONE OF BETTER/WORSE WITH 50/50 PROBABILITY
    h.print_annealing_status();
  }
  
  
  
/*  loop {
    //h.randomize_edge_colors();
    h.randomly_grow_a_clique();
    
    h.find_cliques_from_scratch();
    h.best_anneal_clique_count = h.current_clique_count;

    loops_without_improvement = 0;
    
    while loops_without_improvement < h.graph_size  {
      loops_without_improvement += 1;

      // Try all alternate colors for the current edge
      for _ in 0..(color_count - 1) {
        h.increment_edge_color(edge_index_to_try);
        //h.randomly_grow_a_clique();
        h.find_cliques_from_scratch();
        if h.current_clique_count <= h.best_anneal_clique_count {
          if h.current_clique_count < h.best_global_clique_count || 
          ((h.current_clique_count == h.best_global_clique_count) && h.current_clique_count <= h.upper_bound) {
            
            h.print_header();
            h.best_global_clique_count = h.current_clique_count;
            println!("{}", h);
          }
          if h.current_clique_count < h.best_anneal_clique_count {
            loops_without_improvement = 0;
            h.best_anneal_clique_count = h.current_clique_count
          }
        }
      }
      
      if loops_without_improvement > 0 {
        // No improvement found; revert the color
        h.increment_edge_color(edge_index_to_try);
      }
      
      // Shift to the next edge index
      edge_index_to_try = (edge_index_to_try + 1) % (h.last_random_edge_index + 1);
    }
    
    h.print_annealing_status();
    h.annealing_count += 1;
  }*/
}



struct HyperGraph {
  edges: Vec<Clique>,
  modified_edges: Vec<(usize, u8)>, // (edge_index, color)
  cliques: FxHashMap<u8, Vec<Clique>>,
  members_of_nonmaximal_cliques: FxHashSet<u64>,
  current_clique_count: usize,
  best_anneal_clique_count: usize,
  best_global_clique_count: usize,
  edge_order: u8,
  color_count: u8,
  graph_order: u8,
  graph_size: usize,
  last_random_edge_index: usize,
  set_bit_getter: math::SetBitGetter,
  util_clique: Clique,
  timer: Instant,
  theorem: String,
  upper_bound: usize,
  annealing_count: usize,
}

impl HyperGraph {
  fn new(edge_order: u8, color_count: u8, graph_order: u8) -> HyperGraph {

    let cliques = FxHashMap::<u8, Vec<Clique>>::default();
    
    let graph_size: usize = math::choose(graph_order as usize, edge_order as usize);
    
    let (upper_bound, theorem) = get_upper_bound(graph_order, edge_order, color_count);
  
    HyperGraph {
      edges: Clique::generate_all_cliques(u8::MAX, edge_order, graph_order),
      modified_edges: Vec::<(usize, u8)>::default(),
      cliques,
      members_of_nonmaximal_cliques: FxHashSet::<u64>::default(),
      current_clique_count: 0,
      best_anneal_clique_count: usize::MAX,
      best_global_clique_count: usize::MAX,
      edge_order,
      color_count,
      graph_order,
      graph_size,
      last_random_edge_index: graph_size - 1, 
      set_bit_getter: math::SetBitGetter::new(),
      util_clique: Clique::new(0, 0, 0),
      timer: Instant::now(),
      theorem,
      upper_bound,
      annealing_count: 0,
    }
  }
  
  fn get_key(&self, color: u8, clique_order: u8) -> u8 {
    color * self.graph_order + clique_order
  }
  
  fn rotate_growth_key(&self, key: &mut u8) {
    let color: u8 = (*key / self.graph_order + 1) % self.color_count;
    let mut clique_order: u8 = *key % self.graph_order;
    if color == 0 {
      clique_order += 1;
      if clique_order > self.graph_order {
        clique_order = self.edge_order - 1;
      }
    }
    *key = self.get_key(color, clique_order);
  }

  fn randomize_edge_colors(&mut self) {
    self.annealing_count += 1;
    for i in 0..=self.last_random_edge_index {
      self.edges[i].set_color(fastrand::u8(0..self.color_count));
    }
  }
  
  fn set_edge_colors(&mut self, target_color: u8) {
    for i in 0..=self.last_random_edge_index {
      self.edges[i].set_color(target_color);
    }
  }
  
  fn randomly_grow_a_clique(&mut self) {
	  let mut clique_to_grow: &Clique = &self.util_clique; // A dummy until we find our clique
    let clique_to_grow_indx: usize = fastrand::usize(0..self.current_clique_count);
    let mut vec_min_index: usize = 0;
    let mut vec_max_index: usize;
    let mut clique_index: usize = 0;
    let mut clique_key: u8 = 0;
	  for (cur_key, cliques_vec) in self.cliques.iter() {
      if cliques_vec.is_empty() {
        continue;
      }
      vec_max_index = vec_min_index + cliques_vec.len() - 1;

      if vec_max_index >= clique_to_grow_indx {
        clique_index = clique_to_grow_indx - vec_min_index;
		    clique_to_grow = &cliques_vec[clique_index];
        
        if clique_to_grow.order == self.graph_order {
          println!("WE ONLY REACH THIS IF WE HIT A CLIQUE CONTAINING ALL VERTICES!");
          return;
        }
        clique_key = *cur_key;
		    break;
	    }
      
      vec_min_index = vec_max_index + 1;
    }
    
    // Randomly pick a vertex not in that clique
    let mut vertex: u64 = fastrand::u64(0..self.graph_order as u64);
    while vertex & clique_to_grow.members != 0 {
      vertex = fastrand::u64(0..self.graph_order as u64);
    }
    
    // Recolor all edges whose vertices are all either our chosen
    // vertex or our chosen clique with the clique's colors.
    self.grow_clique(clique_key, clique_index, vertex)
  }
  


  fn grow_clique(&mut self, key: u8, index: usize, vertex: u64) {
    // Recolor all edges whose vertices are all either our chosen
    // vertex or our chosen clique with the clique's colors.
    let color: u8 = self.cliques[&key][index].color;
    let clique_members_with_vertex: u64 = self.cliques[&key][index].members | vertex;
    self.modified_edges.clear();
    for i in 0..=self.last_random_edge_index {
      if self.edges[i].members & clique_members_with_vertex == self.edges[i].members
      && self.edges[i].color != color {
        self.modified_edges.push((i, self.edges[i].color));
        self.edges[i].set_color(color);
      }
    } 
  }
  
  
  fn revert_clique_growth(&mut self) {
    for (i, color) in &self.modified_edges {
      self.edges[*i].color = *color;
    }
  }
  
  fn increment_edge_color(&mut self, edge_index: usize) {
    self.edges[edge_index].increment_color(self.color_count);
  }
  
  fn confirm_cliques_vec_initialized(&mut self, key: u8) {
    self.cliques.entry(key).or_default();
  }
  
  fn add_clique(&mut self, clique: Clique) {
    self.current_clique_count += 1;
    let key = self.get_key(clique.color, clique.order);
    self.confirm_cliques_vec_initialized(key);
    self.cliques.get_mut(&key).expect("Uninitalized Vector").push(clique);
  }
  
  fn add_clique_vec(&mut self, mut clique_vec: Vec<Clique>, key: u8) {
    self.current_clique_count += clique_vec.len();
    self.confirm_cliques_vec_initialized(key);
    self.cliques.get_mut(&key).expect("Uninitalized Vector").append(&mut clique_vec);
  }
  
  fn deactivate_subcliques(&mut self, color: u8, order_of_larger_clique: u8) {
    let key_big = self.get_key(color, order_of_larger_clique);
    if !&self.cliques.contains_key(&key_big) {
      return;
    }
    
    let key_small = self.get_key(color, order_of_larger_clique - 1);
    
    self.members_of_nonmaximal_cliques.clear();
    
    for big_clique in &self.cliques[&key_big] {
      self.set_bit_getter.reset(big_clique.members);
      while self.set_bit_getter.are_there_more_bits() {
        self.members_of_nonmaximal_cliques.insert(big_clique.members ^ self.set_bit_getter.get_bit());
      }
    }


    // Parse all small cliques, deleting any whose members are marked for deletion.
    let mut num_cliques: usize = self.cliques[&key_small].len();
    let mut i = 0;
    while i < num_cliques {
      if self.members_of_nonmaximal_cliques.contains(&self.cliques[&key_small][i].members) {
        self.members_of_nonmaximal_cliques.remove(&self.cliques[&key_small][i].members);
        self.cliques.get_mut(&key_small).expect("Uninitalized Vector").swap_remove(i);
        num_cliques -= 1;
        self.current_clique_count -= 1;
      } else {
        i += 1;
      }
    }
  }

  
  fn find_cliques_from_scratch(&mut self) {
    self.cliques.clear();
    self.current_clique_count = 0;
    
    let mut new_cliques = Vec::<Clique>::new();
    let mut key_small_cliques: u8;

    for color in 0..self.color_count {
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
      self.deactivate_subcliques(color, self.edge_order);
       
      // 'smaller_cliques' refers to all cliques that match the color
      // we're examining, and have (order-1) vertices.
      // The k'th bit of a cliques_index corresponds to the index in cliques[(color, order - 1)]
      
      for order in (self.edge_order + 1)..(self.graph_order + 1) {
        key_small_cliques = self.get_key(color, order - 1);
        
        // Confirm we have enough smaller cliques to build a bigger clique out of
        let smaller_cliques = match self.cliques.get(&key_small_cliques) {
          Some(cliques) => cliques,
          None => break,
        };
        
        if smaller_cliques.len() < order as usize {
          break;
        }

        new_cliques.clear(); // Clear previous cliques

        let expansions = math::get_expansions(&smaller_cliques.iter().map(|c| c.members).collect(), self.graph_order as u64);
        
        for candidate in expansions {
          new_cliques.push(Clique::new(candidate, color, self.graph_order));
        }


        // Move new_cliques into add_clique
        while let Some(clique) = new_cliques.pop() {
          self.add_clique(clique);
        }


        // Having finished finding color-order cliques, mark any nonmaximal smaller cliques as inactive
        self.deactivate_subcliques(color, order);
      }
    }
  }

  fn print_header(&self) {
    println!(
      "\n------------------------------------------\n\n{}.\n\nIMPROVEMENT: {} -> {}. Annealings: {}, Time: {:?}", 
      self.theorem, 
      self.best_global_clique_count, 
      self.current_clique_count,
      self.annealing_count,
      self.timer.elapsed()
    );
  }
  
  fn print_annealing_status(&self) {
    if (self.annealing_count <=           100)                                       || 
       (self.annealing_count <=         1_000 && self.annealing_count %          10 == 0) ||
       (self.annealing_count <=        10_000 && self.annealing_count %         100 == 0) ||
       (self.annealing_count <=       100_000 && self.annealing_count %       1_000 == 0) ||
       (self.annealing_count <=     1_000_000 && self.annealing_count %      10_000 == 0) ||
       (self.annealing_count <=    10_000_000 && self.annealing_count %     100_000 == 0) ||
       (self.annealing_count <=   100_000_000 && self.annealing_count %   1_000_000 == 0) ||
       (self.annealing_count <= 1_000_000_000 && self.annealing_count %  10_000_000 == 0) ||
       (self.annealing_count >  1_000_000_000 && self.annealing_count % 100_000_000 == 0) {
      println!("  annealings: {}, time: {:?}, time per annealing: {:?}", 
        self.annealing_count, 
        self.timer.elapsed(),
        self.timer.elapsed() / self.annealing_count.try_into().unwrap());
    }  
  }
  
  
  fn get_string(&self) -> String {
    let mut vertex_clique_counts: Vec<u8> = vec![0; (self.color_count * self.graph_order) as usize];
    let mut cliques_str = format!("  cargo run --release {} {} {}", self.edge_order, self.color_count, self.graph_order);
    let mut color_clique_count: usize;
    let mut color_clique_str = format!("color clique counts: ");
    let mut ret_str = format!(
      "\nH: A complete {}-uniform hypergraph on {} vertices with {} colors, {} cliques, {} edges\n", 
      self.edge_order, 
      self.graph_order, 
      self.color_count,
      self.current_clique_count,
      self.graph_size);
      
    // Printing edges is too noisy!  
    /*for edge in &self.edges {
      ret_str += &format!("{}\n", edge);
    }*/

    ret_str += "\nMaximal Color-Cliques\n";
    for color in 0..self.color_count {
      color_clique_count = 0;
      ret_str += &format!("  color {}:\n", color);
      let mut first_of_color: bool = true;
      for order in (self.edge_order - 1)..(self.graph_order + 1) {
        let mut order_str = format!("    order {}:\n", order);
        let mut include_order_str: bool = false;
        if !&self.cliques.contains_key(&self.get_key(color, order)) {
          continue;
        }
        for clique in &self.cliques[&self.get_key(color, order)] {
          color_clique_count += 1;
          if order >= self.edge_order {
            if first_of_color {
              cliques_str += " ";
              first_of_color = false;
            } else {
              cliques_str += ",";
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
      color_clique_str += &format!(" ({}: {})", color, color_clique_count);
    }
    
    ret_str += &format!("\n  {} Maximal Color Cliques Found. Vertices on few cliques: (", self.current_clique_count);
    //ret_str = format!("\n  {} Maximal Color Cliques Found. Vertices on few cliques: (", self.current_clique_count);

    
    for color in 0..self.color_count {
      let mut min_count = self.graph_order;
      for i in (color * self.graph_order)..((color + 1) * self.graph_order) {
        if vertex_clique_counts[i as usize] < min_count {
          min_count = vertex_clique_counts[i as usize];
        }
      }
      ret_str += &format!("{}", min_count);
      if color != self.color_count - 1 {
        ret_str += ", ";
      }
    }
    ret_str += ")\n";
    ret_str += &format!("  {}\n\n",color_clique_str);
    ret_str += &cliques_str;
    ret_str += "\n";
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
}

impl Clique {
  fn new(members: u64, color: u8, graph_order: u8) -> Clique {
    Clique {
      members,
      color,
      order: members.count_ones() as u8,
      graph_order,
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
  
  fn increment_color(&mut self, color_count: u8) {
    self.color = (self.color + 1) % color_count;
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
  
  if r == 2 && t_is_ppo_plus_one {
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
    } else if n_mod_5 == 2 || n_mod_5 == 3 {
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

fn lemma9_n1000(n: usize, c0: usize, c1: usize, d0: usize, d1: usize) -> usize {
  // ci = color i cliques; di = number of color-i cliques for vertex on fewest of them
  if d0 > d1 {
    return lemma9_n1000(n, c1, c0, d1, d0);
  }
  let mut cliques: usize = c0 + c1;
  if d0 < d1 {
    cliques += (d1-1)*(d1)/2;
    cliques -= (d0-1)*(d0)/2;
  }
  let pairwise_steps: usize = (1000 - (n + d1 - d0))/2;
  cliques += pairwise_steps * (pairwise_steps + 1); // we don't divide by 2 b/c we do this twice
  cliques -= (d1 - 1) * (d1);

  // If we only have 999 vertices, add one more
  if (n + d1 - d0) % 2 == 1 {
    cliques += pairwise_steps + d1;
  }
  return cliques;
}
