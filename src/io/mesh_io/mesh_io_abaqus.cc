/**
 * @file   mesh_io_abaqus.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Fri Dec 11 2015
 *
 * @brief  read a mesh from an abaqus input file
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
// std library header files
#include <fstream>

// akantu header files
#include "mesh_io_abaqus.hh"
#include "mesh.hh"
#include "mesh_utils.hh"

#include "element_group.hh"
#include "node_group.hh"


#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
#elif (defined(__GNUC__) || defined(__GNUG__))
#  define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#  if GCC_VERSION > 40600
#    pragma GCC diagnostic push
#  endif
#  pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

/* -------------------------------------------------------------------------- */
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
MeshIOAbaqus::MeshIOAbaqus() {}

/* -------------------------------------------------------------------------- */
MeshIOAbaqus::~MeshIOAbaqus() {}

/* -------------------------------------------------------------------------- */
namespace spirit = boost::spirit;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace lbs = boost::spirit::qi::labels;
namespace phx = boost::phoenix;

/* -------------------------------------------------------------------------- */
void element_read(Mesh & mesh, const ElementType & type, UInt id, const std::vector<Int> & conn,
                  const std::map<UInt, UInt> & nodes_mapping,
                  std::map<UInt, Element> & elements_mapping) {
  Vector<UInt> tmp_conn(Mesh::getNbNodesPerElement(type));

  AKANTU_DEBUG_ASSERT(conn.size() == tmp_conn.size(),
		      "The nodes in the Abaqus file have too many coordinates"
		      << " for the mesh you try to fill.");

  mesh.addConnectivityType(type);
  Array<UInt> & connectivity = mesh.getConnectivity(type);

  UInt i = 0;
  for (std::vector<Int>::const_iterator it = conn.begin(); it != conn.end(); ++it) {
    std::map<UInt, UInt>::const_iterator nit = nodes_mapping.find(*it);
    AKANTU_DEBUG_ASSERT(nit != nodes_mapping.end(),
			"There is an unknown node in the connectivity.");
    tmp_conn[i++] = nit->second;
  }
  Element el(type, connectivity.getSize());
  elements_mapping[id] = el;
  connectivity.push_back(tmp_conn);
}


void node_read(Mesh & mesh, UInt id, const std::vector<Real> & pos,
	       std::map<UInt, UInt> & nodes_mapping) {
  Vector<Real> tmp_pos(mesh.getSpatialDimension());
  UInt i = 0;
  for (std::vector<Real>::const_iterator it = pos.begin();
       it != pos.end() || i < mesh.getSpatialDimension(); ++it)
    tmp_pos[i++] = *it;

  nodes_mapping[id] = mesh.getNbNodes();
  mesh.getNodes().push_back(tmp_pos);
}


/* ------------------------------------------------------------------------ */
void add_element_to_group(ElementGroup * el_grp, UInt element,
			  const std::map<UInt, Element> & elements_mapping) {
  std::map<UInt, Element>::const_iterator eit = elements_mapping.find(element);
  AKANTU_DEBUG_ASSERT(eit != elements_mapping.end(),
		      "There is an unknown element ("
		      << element << ") in the in the ELSET "
		      << el_grp->getName() << ".");

  el_grp->add(eit->second, true, false);
}

ElementGroup * element_group_create(Mesh & mesh, const ID & name) {
  Mesh::element_group_iterator eg_it = mesh.element_group_find(name);
  if (eg_it != mesh.element_group_end()) {
    return eg_it->second;
  } else {
    return &mesh.createElementGroup(name, _all_dimensions);
  }
}

NodeGroup * node_group_create(Mesh & mesh, const ID & name) {
  Mesh::node_group_iterator ng_it = mesh.node_group_find(name);
  if (ng_it != mesh.node_group_end()) {
    return ng_it->second;
  } else {
    return &mesh.createNodeGroup(name, mesh.getSpatialDimension());
  }
}

void add_node_to_group(NodeGroup * node_grp, UInt node,
		       const std::map<UInt, UInt> & nodes_mapping) {
  std::map<UInt, UInt>::const_iterator nit = nodes_mapping.find(node);

  AKANTU_DEBUG_ASSERT(nit != nodes_mapping.end(),
		      "There is an unknown node in the in the NSET "
		      << node_grp->getName() << ".");
  
  node_grp->add(nit->second, false);
}

void optimize_group(NodeGroup * grp) { grp->optimize(); }
void optimize_element_group(ElementGroup * grp) { grp->optimize(); }

/* -------------------------------------------------------------------------- */
template <class Iterator> struct AbaqusSkipper : qi::grammar<Iterator> {
  AbaqusSkipper() : AbaqusSkipper::base_type(skip, "abaqus_skipper") {
    /* clang-format off */
    skip
      =   (ascii::space - spirit::eol)
      |   "**" >> *(qi::char_ - spirit::eol) >> spirit::eol
      ;
    /* clang-format on */
  }
  qi::rule<Iterator> skip;
};

/* -------------------------------------------------------------------------- */
template <class Iterator, typename Skipper = AbaqusSkipper<Iterator> >
struct AbaqusMeshGrammar : qi::grammar<Iterator, void(), Skipper> {
public:
  AbaqusMeshGrammar(Mesh & mesh)
      : AbaqusMeshGrammar::base_type(start, "abaqus_mesh_reader"), mesh(mesh) {

    /* clang-format off */
    start
      =  *(
             (qi::char_('*')
                >  (   (qi::no_case[ qi::lit("node output")    ] > any_section)
                   |   (qi::no_case[ qi::lit("element output") ] > any_section)
                   |   (qi::no_case[ qi::lit("node")           ] > nodes)
                   |   (qi::no_case[ qi::lit("element")        ] > elements)
                   |   (qi::no_case[ qi::lit("heading")        ] > header)
                   |   (qi::no_case[ qi::lit("elset")          ] > elements_set)
                   |   (qi::no_case[ qi::lit("nset")           ] > nodes_set)
                   |   (qi::no_case[ qi::lit("material")       ] > material)
                   |   (keyword > any_section)
                   )
             )
          |  spirit::eol
          )
      ;

    header
      =   spirit::eol
      >   *any_line
      ;

    nodes
      =   *(qi::char_(',') >> option)
           >> spirit::eol
           >> *( (qi::int_
                  > node_position) [ phx::bind(&node_read,
						 phx::ref(mesh),
						 lbs::_1,
						 lbs::_2,
						 phx::ref(abaqus_nodes_to_akantu)) ]
                  >> spirit::eol
               )
      ;

    elements
      =   (
             (  qi::char_(',') >> qi::no_case[qi::lit("type")] >> '='
                >> abaqus_element_type  [ lbs::_a = lbs::_1 ]
             )
          ^  *(qi::char_(',') >> option)
          )
          >> spirit::eol
          >> *(  (qi::int_
                  > connectivity) [ phx::bind(&element_read,
						phx::ref(mesh),
						lbs::_a,
						lbs::_1,
						lbs::_2,
						phx::cref(abaqus_nodes_to_akantu),
						phx::ref(abaqus_elements_to_akantu)) ]
                 >> spirit::eol
               )
      ;

    elements_set
      =   (
             (
                (  qi::char_(',') >> qi::no_case[ qi::lit("elset") ] >> '='
                   >> value [ lbs::_a = phx::bind<ElementGroup *>(&element_group_create,
								    phx::ref(mesh),
								    lbs::_1) ]
                )
             ^  *(qi::char_(',') >> option)
             )
             >> spirit::eol
             >> qi::skip
                  (qi::char_(',') | qi::space)
	     [ +(qi::int_ [ phx::bind(&add_element_to_group,
					lbs::_a,
					lbs::_1,
					phx::cref(abaqus_elements_to_akantu)
					)
                               ]
                     )
                 ]
	   ) [ phx::bind(&optimize_element_group, lbs::_a) ]
      ;

    nodes_set
      =   (
             (
                (  qi::char_(',')
                   >> qi::no_case[ qi::lit("nset") ] >> '='
                   >> value [ lbs::_a = phx::bind<NodeGroup *>(&node_group_create, phx::ref(mesh), lbs::_1) ]
                )
             ^  *(qi::char_(',') >> option)
             )
             >> spirit::eol
             >> qi::skip
                 (qi::char_(',') | qi::space)
	         [ +(qi::int_ [ phx::bind(&add_node_to_group,
					    lbs::_a,
					    lbs::_1,
					    phx::cref(abaqus_nodes_to_akantu)
					   )
                              ]
                    )
                 ]
           ) [ phx::bind(&optimize_group, lbs::_a) ]
      ;

    material
      =   (
             (  qi::char_(',') >> qi::no_case[ qi::lit("name") ] >> '='
                >> value  [ phx::push_back(phx::ref(material_names), lbs::_1) ]
             )
          ^  *(qi::char_(',') >> option)
          ) >> spirit::eol;
      ;

    node_position
      =   +(qi::char_(',') > real [ phx::push_back(lbs::_val, lbs::_1) ])
      ;

    connectivity
      = +(qi::char_(',') > qi::int_ [ phx::push_back(lbs::_val, lbs::_1) ])
      ;


    any_section
      =   *(qi::char_(',') >> option) > spirit::eol
      >   *any_line
      ;

    any_line
      =   *(qi::char_ - spirit::eol - qi::char_('*')) >> spirit::eol
      ;

    keyword
      =   qi::lexeme[ +(qi::char_ - (qi::char_('*') | spirit::eol)) ]
      ;

    option
      = key > -( '=' >> value );

    key
      =   qi::char_("a-zA-Z_") >> *(qi::char_("a-zA-Z_0-9") | qi::char_('-'))
      ;

    value
      =   key.alias()
      ;

    BOOST_SPIRIT_DEBUG_NODE(start);

    abaqus_element_type.add
#if defined(AKANTU_STRUCTURAL_MECHANICS)
      ("BE21"  , _bernoulli_beam_2)
      ("BE31"  , _bernoulli_beam_3)
#endif
      ("T3D2"  , _segment_2) // Gmsh generates this elements
      ("T3D3"  , _segment_3) // Gmsh generates this elements
      ("CPE3"  , _triangle_3)
      ("CPS3"  , _triangle_3)
      ("DC2D3" , _triangle_3)
      ("CPE6"  , _triangle_6)
      ("CPS6"  , _triangle_6)
      ("DC2D6" , _triangle_6)
      ("CPE4"  , _quadrangle_4)
      ("CPS4"  , _quadrangle_4)
      ("DC2D4" , _quadrangle_4)
      ("CPE8"  , _quadrangle_8)
      ("CPS8"  , _quadrangle_8)
      ("DC2D8" , _quadrangle_8)
      ("C3D4"  , _tetrahedron_4)
      ("DC3D4" , _tetrahedron_4)
      ("C3D8"  , _hexahedron_8)
      ("C3D8R" , _hexahedron_8)
      ("DC3D8" , _hexahedron_8)
      ("C3D10" , _tetrahedron_10)
      ("DC3D10", _tetrahedron_10);

#if !defined(AKANTU_NDEBUG) && defined(AKANTU_CORE_CXX_11)
    qi::on_error<qi::fail>(start, error_handler(lbs::_4, lbs::_3, lbs::_2));
#endif
    
    start              .name("abaqus-start-rule");
    connectivity       .name("abaqus-connectivity");
    node_position      .name("abaqus-nodes-position");
    nodes              .name("abaqus-nodes");
    any_section        .name("abaqus-any_section");
    header             .name("abaqus-header");
    material           .name("abaqus-material");
    elements           .name("abaqus-elements");
    elements_set       .name("abaqus-elements-set");
    nodes_set          .name("abaqus-nodes-set");
    key                .name("abaqus-key");
    value              .name("abaqus-value");
    option             .name("abaqus-option");
    keyword            .name("abaqus-keyword");
    any_line           .name("abaqus-any-line");
    abaqus_element_type.name("abaqus-element-type");

    /* clang-format on */
  }

public:
  AKANTU_GET_MACRO(MaterialNames, material_names,
                   const std::vector<std::string> &);

  /* ------------------------------------------------------------------------ */
  /* Rules                                                                    */
  /* ------------------------------------------------------------------------ */
private:
  qi::rule<Iterator, void(), Skipper> start;
  qi::rule<Iterator, std::vector<int>(), Skipper> connectivity;
  qi::rule<Iterator, std::vector<Real>(), Skipper> node_position;
  qi::rule<Iterator, void(), Skipper> nodes, any_section, header, material;
  qi::rule<Iterator, void(), qi::locals<ElementType>, Skipper> elements;
  qi::rule<Iterator, void(), qi::locals<ElementGroup *>, Skipper> elements_set;
  qi::rule<Iterator, void(), qi::locals<NodeGroup *>, Skipper> nodes_set;

  qi::rule<Iterator, std::string(), Skipper> key, value, option, keyword,
      any_line;

  qi::real_parser<Real, qi::real_policies<Real> > real;

  qi::symbols<char, ElementType> abaqus_element_type;

  /* ------------------------------------------------------------------------ */
  /* Mambers                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  /// reference to the mesh to read
  Mesh & mesh;

  /// correspondance between the numbering of nodes in the abaqus file and in
  /// the akantu mesh
  std::map<UInt, UInt> abaqus_nodes_to_akantu;

  /// correspondance between the element number in the abaqus file and the
  /// Element in the akantu mesh
  std::map<UInt, Element> abaqus_elements_to_akantu;

  /// list of the material names
  std::vector<std::string> material_names;
};
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void MeshIOAbaqus::read(const std::string & filename, Mesh & mesh) {
  namespace spirit = boost::spirit;
  namespace qi = boost::spirit::qi;
  namespace lbs = boost::spirit::qi::labels;
  namespace ascii = boost::spirit::ascii;
  namespace phx = boost::phoenix;

  std::ifstream infile;
  infile.open(filename.c_str());

  if (!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }

  std::string storage;             // We will read the contents here.
  infile.unsetf(std::ios::skipws); // No white space skipping!
  std::copy(std::istream_iterator<char>(infile), std::istream_iterator<char>(),
            std::back_inserter(storage));

  typedef std::string::const_iterator iterator_t;
  typedef AbaqusSkipper<iterator_t> skipper;
  typedef AbaqusMeshGrammar<iterator_t, skipper> grammar;

  grammar g(mesh);
  skipper ws;

  iterator_t iter = storage.begin();
  iterator_t end = storage.end();

  qi::phrase_parse(iter, end, g, ws);

  std::vector<std::string>::const_iterator mnit = g.getMaterialNames().begin();
  std::vector<std::string>::const_iterator mnend = g.getMaterialNames().end();

  MeshAccessor mesh_accessor(mesh);

  for (; mnit != mnend; ++mnit) {
    Mesh::element_group_iterator eg_it = mesh.element_group_find(*mnit);
    ElementGroup & eg = *eg_it->second;
    if (eg_it != mesh.element_group_end()) {
      ElementGroup::type_iterator tit = eg.firstType();
      ElementGroup::type_iterator tend = eg.lastType();

      for (; tit != tend; ++tit) {
        Array<std::string> & abaqus_material =
            mesh_accessor.getData<std::string>("abaqus_material", *tit);

        ElementGroup::const_element_iterator eit = eg.element_begin(*tit);
        ElementGroup::const_element_iterator eend = eg.element_end(*tit);
        for (; eit != eend; ++eit) {
          abaqus_material(*eit) = *mnit;
        }
      }
    }
  }

  mesh_accessor.setNbGlobalNodes(mesh.getNodes().getSize());
  MeshUtils::fillElementToSubElementsData(mesh);
}

} // akantu
