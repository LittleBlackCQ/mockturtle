#include <lorina/aiger.hpp>
#include <mockturtle/mockturtle.hpp>
#include <mockturtle/algorithms/extract_adders.hpp>
#include <mockturtle/io/write_dot.hpp>

int main(int argc, char* argv[])
{
  if(argc != 2 && argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input_aig_file> [-o dot_file_name]" << std::endl;
    return 1;
  }
  std::string const filename = argv[1];

  mockturtle::aig_network aig;
  auto const result = lorina::read_aiger( filename, mockturtle::aiger_reader( aig ) );
  
  if (result != lorina::return_code::success) {
    std::cerr << "Error: Failed to read AIG file: " << filename << std::endl;
    return 1;
  }

  mockturtle::extract_adders_params params;
  params.verbose = true;
  auto block_network = extract_adders(aig, params);

  if (argc == 4) {
    std::string option( argv[2] );
    std::string dot_path ( argv[3] );
    if (option != "-o") {
      std::cerr << "Usage: " << argv[0] << " <input_aig_file> [-o dot_file_name]" << std::endl;
      return 1;
    }
    mockturtle::write_dot(block_network, dot_path, mockturtle::adders_dot_drawer<mockturtle::block_network>());
  }
  return 0;
}
