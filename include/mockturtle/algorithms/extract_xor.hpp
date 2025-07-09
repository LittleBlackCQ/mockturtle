/*!
  \file extract_xor.hpp
  \brief Extracts XOR gates in the network

  \author Generated based on extract_adders.hpp
*/

#include <algorithm>
#include <array>
#include <vector>

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <parallel_hashmap/phmap.h>

#include "../networks/block.hpp"
#include "../networks/storage.hpp"
#include "../utils/node_map.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/choice_view.hpp"
#include "cut_enumeration.hpp"

namespace mockturtle
{

struct extract_xor_params
{
  extract_xor_params()
  {
    cut_enumeration_ps.cut_limit = 49;
    cut_enumeration_ps.minimize_truth_table = false;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut limit is 49. By default,
   * truth table minimization is not performed.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Be verbose */
  bool verbose{ false };
};

struct extract_xor_stats
{
  /*! \brief Computed cuts. */
  uint32_t cuts_total{ 0 };

  /*! \brief XOR2 gates count. */
  uint32_t xor2{ 0 };

  /*! \brief Mapped XOR2 gates. */
  uint32_t mapped_xor2{ 0 };

  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] Cuts = {}\t XOR2 = {}\n", cuts_total, xor2 );
    std::cout << fmt::format( "[i] Mapped XOR2 = {}\n", mapped_xor2 );
    std::cout << fmt::format( "[i] Total runtime = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

struct double_hash
{
  uint64_t operator()( const std::array<uint32_t, 2>& p ) const
  {
    uint64_t seed = hash_block( p[0] );
    hash_combine( seed, hash_block( p[1] ) );
    return seed;
  }
};

struct cut_enumeration_xor_cut
{
  /* no additional data needed for XOR extraction */
};

template<class Ntk>
class extract_xor_impl
{
public:
  using network_cuts_t = fast_network_cuts<Ntk, 2, true, cut_enumeration_xor_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using leaves_hash_t = phmap::flat_hash_map<std::array<uint32_t, 2>, std::vector<uint64_t>, double_hash>;
  using block_map = node_map<signal<block_network>, Ntk>;

public:
  explicit extract_xor_impl( Ntk& ntk, extract_xor_params const& ps, extract_xor_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st ),
        cuts( fast_cut_enumeration<Ntk, 2, true, cut_enumeration_xor_cut>( ntk, ps.cut_enumeration_ps ) ),
        cuts_classes(),
        xor_candidates(),
        node_match( ntk.size(), UINT32_MAX )
  {
    cuts_classes.reserve( 1000 );
  }

  block_network run()
  {
    stopwatch t( st.time_total );

    auto [res, old2new] = initialize_map_network();
    create_classes();
    match_xor_gates();
    map();
    topo_sort();
    finalize( res, old2new );

    return res;
  }

private:
  void create_classes()
  {
    st.cuts_total = cuts.total_cuts();

    ntk.foreach_gate( [&]( auto const& n ) {
      uint32_t cut_index = 0;
      for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
      {
        if ( cut->size() != 2 )
        {
          ++cut_index;
          continue;
        }

        kitty::static_truth_table<2> tt = cuts.truth_table( *cut );

        /* check for xor2 */
        bool is_xor = false;
        for ( uint32_t func : xor2func )
        {
          if ( tt._bits == func )
          {
            ++st.xor2;
            is_xor = true;
            break;
          }
        }

        if ( !is_xor )
        {
          ++cut_index;
          continue;
        }

        uint64_t data = ( static_cast<uint64_t>( ntk.node_to_index( n ) ) << 16 ) | cut_index;
        std::array<uint32_t, 2> leaves = { 0, 0 };
        uint32_t i = 0;
        for ( auto l : *cut )
          leaves[i++] = l;

        /* sort leaves to ensure consistent ordering */
        if ( leaves[0] > leaves[1] )
          std::swap( leaves[0], leaves[1] );

        /* add to hash table */
        auto& v = cuts_classes[leaves];
        v.push_back( data );

        ++cut_index;
      }
    } );
  }

  void match_xor_gates()
  {
    xor_candidates.reserve( cuts_classes.size() );

    for ( auto& it : cuts_classes )
    {
      /* For XOR gates, we can have multiple equivalent implementations */
      for ( uint32_t i = 0; i < it.second.size(); ++i )
      {
        uint64_t data = it.second[i];
        xor_candidates.push_back( data );
      }
    }
  }

  void map()
  {
    selected.reserve( xor_candidates.size() );

    ntk.incr_trav_id();

    for ( uint32_t i = 0; i < xor_candidates.size(); ++i )
    {
      uint64_t data = xor_candidates[i];
      uint32_t node_index = data >> 16;

      /* check if node is already mapped */
      if ( node_match[node_index] != UINT32_MAX )
        continue;

      selected.push_back( i );
      node_match[node_index] = i;
      ++st.mapped_xor2;
    }
  }

  void topo_sort()
  {
    topo_order.reserve( ntk.size() );

    ntk.incr_trav_id();
    ntk.incr_trav_id();

    /* add constants and CIs */
    const auto c0 = ntk.get_node( ntk.get_constant( false ) );
    ntk.set_visited( c0, ntk.trav_id() );

    if ( const auto c1 = ntk.get_node( ntk.get_constant( true ) ); ntk.visited( c1 ) != ntk.trav_id() )
    {
      ntk.set_visited( c1, ntk.trav_id() );
    }

    ntk.foreach_ci( [&]( auto const& n ) {
      if ( ntk.visited( n ) != ntk.trav_id() )
      {
        ntk.set_visited( n, ntk.trav_id() );
      }
    } );

    /* sort topologically */
    ntk.foreach_co( [&]( auto const& f ) {
      if ( ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() )
        return;
      topo_sort_rec( ntk.get_node( f ) );
    } );
  }
  
  bool is_xor_candidate( node<Ntk> const& n )
  {
    for (auto& data : xor_candidates)
    {
      uint32_t node_index = data >> 16;
      if (ntk.node_to_index(n) == node_index)
        return true;
    }
    return false;
  }

  void topo_sort_rec( node<Ntk> const& n )
  {
    /* is permanently marked? */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    if (is_xor_candidate(n))
    {
      ntk.set_visited( n, ntk.trav_id() - 1 );
      uint32_t node_index = ntk.node_to_index(n);
      auto& data = xor_candidates[node_match[node_index]];
      cut_t const& cut = cuts.cuts( data >> 16 )[data & UINT16_MAX];

      for ( auto l : cut )
      {
        topo_sort_rec( ntk.index_to_node(l) );
      }
      ntk.set_visited( n, ntk.trav_id() );
      topo_order.push_back( n );
    }
    else
    { 
      /* ensure that the node is not visited or temporarily marked */
      assert( ntk.visited( n ) != ntk.trav_id() );
      assert( ntk.visited( n ) != ntk.trav_id() - 1 );

      /* mark node temporarily */
      ntk.set_visited( n, ntk.trav_id() - 1 );

      /* mark cut leaves */
      ntk.foreach_fanin( n, [&]( auto const& f ) {
      topo_sort_rec( ntk.get_node( f ) );
      } );

      /* ensure that the node is not visited */
      assert( ntk.visited( n ) != ntk.trav_id() );

      /* mark node n permanently */
      ntk.set_visited( n, ntk.trav_id() );

      /* visit node */
      topo_order.push_back( n );
    }
  }

  std::pair<block_network, block_map> initialize_map_network()
  {
    block_network dest;
    block_map old2new( ntk );

    old2new[ntk.get_node( ntk.get_constant( false ) )] = dest.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
      old2new[ntk.get_node( ntk.get_constant( true ) )] = dest.get_constant( true );

    ntk.foreach_ci( [&]( auto const& n ) {
      old2new[n] = dest.create_pi();
    } );
    return { dest, old2new };
  }

  void finalize( block_network& res, block_map& old2new )
  {
    for ( auto const& n : topo_order )
    {
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        continue;

      uint32_t node_index = ntk.node_to_index( n );
      if ( node_match[node_index] != UINT32_MAX )
      {
        /* This is a mapped XOR gate */
        finalize_xor_gate( res, old2new, n );
      }
      else
      {
        /* This is a regular gate */
        finalize_simple_gate( res, old2new, n );
      }
    }

    /* Create POs */
    ntk.foreach_co( [&]( auto const& f ) {
      res.create_po( ntk.is_complemented( f ) ? !old2new[f] : old2new[f] );
    } );
  }

  inline void finalize_simple_gate( block_network& res, block_map& old2new, node<Ntk> const& n )
  {
    std::vector<signal<block_network>> children;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto s = old2new[f] ^ ntk.is_complemented( f );
      children.push_back( s );
    } );

    /* Create a generic gate with the same function */
    old2new[n] = res.create_node( children, ntk.node_function( n ) );
  }

  inline void finalize_xor_gate( block_network& res, block_map& old2new, node<Ntk> const& n )
  {
    uint32_t node_index = ntk.node_to_index( n );
    uint32_t match_index = node_match[node_index];
    uint64_t data = xor_candidates[match_index];
    uint32_t cut_index = data & UINT16_MAX;

    /* Get the original cut */
    cut_t const& cut = cuts.cuts( node_index )[cut_index];
    kitty::static_truth_table<2> tt = cuts.truth_table( cut );

    /* Find the correct XOR function form and polarity */
    uint32_t func_index = 0;
    for ( uint32_t func : xor2func )
    {
      if ( tt._bits == func )
        break;
      ++func_index;
    }

    /* Create the XOR gate */
    std::array<signal<block_network>, 2> children;
    uint32_t i = 0;
    for ( auto l : cut )
    {
      signal<block_network> f = old2new[ntk.index_to_node( l )];
      /* Apply appropriate polarity based on the XOR function */
      bool polarity = ( func_index == 1 ) && ( i == 0 ); /* XNOR case */
      children[i] = f ^ polarity;
      ++i;
    }

    /* Create XOR2 gate */
    signal<block_network> xor_gate = res.create_xor( children[0], children[1] );
    /* Apply output polarity if needed */
    old2new[n] = xor_gate ^ ( func_index == 1 );
  }

private:
  Ntk& ntk;
  extract_xor_params const& ps;
  extract_xor_stats& st;

  network_cuts_t cuts;
  leaves_hash_t cuts_classes;
  std::vector<uint64_t> xor_candidates;
  std::vector<uint32_t> selected;
  std::vector<uint32_t> node_match;

  std::vector<node<Ntk>> topo_order;
  /* XOR2 truth table functions: XOR and XNOR */
  const std::array<uint32_t, 2> xor2func = { 0x6, 0x9 };  // XOR: 0110, XNOR: 1001
};

} /* namespace detail */

/*! \brief XOR extraction.
 *
 * This function extracts 2-input XOR gates from a network.
 * It returns a `block_network` with extracted XOR gates replaced
 * by dedicated XOR blocks.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_co`
 * - `foreach_node`
 * - `foreach_gate`
 * - `foreach_ci`
 * - `foreach_fanin`
 * - `node_function`
 * - `is_complemented`
 * - `get_constant`
 * - `create_node`
 * - `create_pi`
 * - `create_po`
 * - `create_xor`
 *
 * \param ntk Network
 * \param ps Parameters
 * \param pst Stats
 *
 */
template<class Ntk>
block_network extract_xor( Ntk& ntk, extract_xor_params const& ps = {}, extract_xor_stats* pst = {} )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_foreach_co_v<Ntk>, "Ntk does not implement the foreach_co method" );

  extract_xor_stats st;

  detail::extract_xor_impl p( ntk, ps, st );
  block_network res = p.run();

  if ( ps.verbose )
    st.report();

  if ( pst )
    *pst = st;

  return res;
}

}; /* namespace mockturtle */
