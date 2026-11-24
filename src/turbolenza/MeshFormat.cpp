// The "turbolenza.mesh" JSON format, version 1 (see MeshFormat.hpp and docs/turbolenza/mesh_format.md).
#include "turbolenza/MeshFormat.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <format>
#include <fstream>
#include <initializer_list>
#include <set>
#include <sstream>

namespace BladeGen::Turbolenza {

namespace {

using nlohmann::json;

// JSON string escaping for names and descriptions.
std::string quoted(const std::string& s)
{
    std::string out = "\"";
    for (char c : s) {
        switch (c) {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n"; break;
        case '\t': out += "\\t"; break;
        case '\r': out += "\\r"; break;
        default:
            if (static_cast<unsigned char>(c) < 0x20) out += std::format("\\u{:04x}", static_cast<unsigned>(c));
            else out += c;
        }
    }
    return out + "\"";
}

std::string number(double v) { return std::format("{:.17g}", v); }

std::string array(const std::vector<double>& v)
{
    std::string out = "[";
    for (std::size_t k = 0; k < v.size(); ++k) {
        if (k) out += k % 8 == 0 ? ",\n        " : ", ";
        out += number(v[k]);
    }
    return out + "]";
}

[[noreturn]] void fail(const std::string& what) { throw MeshFormatError("turbolenza.mesh: " + what); }

void only_keys(const json& j, std::initializer_list<const char*> keys, const std::string& where)
{
    if (!j.is_object()) fail(where + " must be an object");
    std::set<std::string> allowed(keys.begin(), keys.end());
    for (auto it = j.begin(); it != j.end(); ++it)
        if (!allowed.count(it.key())) fail(std::format("{} has the unknown key '{}'", where, it.key()));
}

const json& need(const json& j, const char* key, const std::string& where)
{
    if (!j.contains(key)) fail(std::format("{} has no '{}'", where, key));
    return j.at(key);
}

int integer(const json& j, const char* key, const std::string& where)
{
    const json& v = need(j, key, where);
    if (!v.is_number_integer()) fail(std::format("{}: '{}' must be an integer", where, key));
    return v.get<int>();
}

double real(const json& v, const std::string& where)
{
    if (!v.is_number()) fail(where + " must be a number");
    const double x = v.get<double>();
    if (!std::isfinite(x)) fail(where + " is not finite");
    return x;
}

std::string text(const json& j, const char* key, const std::string& where)
{
    const json& v = need(j, key, where);
    if (!v.is_string()) fail(std::format("{}: '{}' must be a string", where, key));
    return v.get<std::string>();
}

std::vector<double> reals(const json& j, const char* key, const std::string& where)
{
    const json& v = need(j, key, where);
    if (!v.is_array()) fail(std::format("{}: '{}' must be an array of numbers", where, key));
    std::vector<double> out;
    out.reserve(v.size());
    for (std::size_t k = 0; k < v.size(); ++k) out.push_back(real(v[k], std::format("{} '{}'[{}]", where, key, k)));
    return out;
}

BlockFace face_of(const std::string& s, const std::string& where)
{
    for (BlockFace f : {BlockFace::IMin, BlockFace::IMax, BlockFace::JMin, BlockFace::JMax})
        if (s == block_face_name(f)) return f;
    fail(std::format("{}: face '{}' is not imin, imax, jmin or jmax", where, s));
}

Surface surface_of(const std::string& s, const std::string& where)
{
    if (s == "suction") return Surface::Suction;
    if (s == "pressure") return Surface::Pressure;
    fail(std::format("{}: surface '{}' is not suction or pressure", where, s));
}

}  // namespace

std::string write_mesh_json(const GeneralMesh& mesh)
{
    validate_mesh(mesh);
    for (const FaceSegment& g : mesh.segments)
        if (g.kind == SegmentKind::Dirichlet)
            fail("a Dirichlet segment (verification only) cannot be written: it is not part of the format");
    std::ostringstream o;
    o << "{\n";
    o << "  \"format\": " << quoted(kMeshFormatName) << ",\n";
    o << "  \"version\": " << kMeshFormatVersion << ",\n";
    o << "  \"units\": {\"m_prime\": \"1\", \"theta\": \"rad\", \"m\": \"m\", \"r\": \"m\", \"b\": \"1\"},\n";
    o << "  \"description\": " << quoted(mesh.description) << ",\n";
    o << "  \"pitch\": " << number(mesh.pitch) << ",\n";
    o << "  \"blocks\": [\n";
    for (std::size_t k = 0; k < mesh.blocks.size(); ++k) {
        const MeshBlock& B = mesh.blocks[k];
        o << "    {\"name\": " << quoted(B.name) << ", \"ni\": " << B.ni << ", \"nj\": " << B.nj << ",\n";
        o << "      \"m_prime\": " << array(B.x) << ",\n";
        o << "      \"theta\": " << array(B.y) << ",\n";
        o << "      \"m\": " << array(B.m) << ",\n";
        o << "      \"r\": " << array(B.r) << ",\n";
        o << "      \"b\": " << array(B.b) << "}" << (k + 1 < mesh.blocks.size() ? "," : "") << "\n";
    }
    o << "  ],\n  \"segments\": [\n";
    for (std::size_t k = 0; k < mesh.segments.size(); ++k) {
        const FaceSegment& g = mesh.segments[k];
        o << "    {\"block\": " << g.block << ", \"face\": " << quoted(block_face_name(g.face)) << ", \"range\": [" << g.begin
          << ", " << g.end << "], \"kind\": " << quoted(segment_kind_name(g.kind));
        if (g.kind == SegmentKind::Blade)
            o << ", \"blade\": " << g.blade << ", \"surface\": " << quoted(mesh_surface_name(g.surface));
        o << "}" << (k + 1 < mesh.segments.size() ? "," : "") << "\n";
    }
    o << "  ],\n  \"connections\": [\n";
    for (std::size_t k = 0; k < mesh.connections.size(); ++k) {
        const Connection& c = mesh.connections[k];
        const bool periodic = mesh.segments[static_cast<std::size_t>(c.first)].kind == SegmentKind::Periodic;
        o << "    {\"kind\": " << quoted(periodic ? "periodic" : "interface") << ", \"segments\": [" << c.first << ", "
          << c.second << "], \"reversed\": " << (c.reversed ? "true" : "false");
        if (periodic) o << ", \"image\": " << c.image;
        o << "}" << (k + 1 < mesh.connections.size() ? "," : "") << "\n";
    }
    o << "  ],\n  \"tags\": [\n";
    for (std::size_t k = 0; k < mesh.tags.size(); ++k) {
        const EdgeTag& t = mesh.tags[k];
        o << "    {\"edge\": " << quoted(edge_kind_name(t.edge)) << ", \"blade\": " << t.blade
          << ", \"surface\": " << quoted(mesh_surface_name(t.surface)) << ", \"block\": " << t.block << ", \"node\": ["
          << t.i << ", " << t.j << "]}" << (k + 1 < mesh.tags.size() ? "," : "") << "\n";
    }
    o << "  ]\n}\n";
    return o.str();
}

GeneralMesh read_mesh_json(const std::string& source)
{
    json root;
    try {
        root = json::parse(source);
    } catch (const json::exception& e) {
        fail(std::string("not JSON: ") + e.what());
    }
    only_keys(root, {"format", "version", "units", "description", "pitch", "blocks", "segments", "connections", "tags"},
              "the file");
    if (text(root, "format", "the file") != kMeshFormatName)
        fail(std::format("the format is '{}', not '{}'", root.at("format").get<std::string>(), kMeshFormatName));
    const int version = integer(root, "version", "the file");
    if (version != kMeshFormatVersion)
        fail(std::format("version {} is not supported (this reader reads version {})", version, kMeshFormatVersion));
    const json& units = need(root, "units", "the file");
    only_keys(units, {"m_prime", "theta", "m", "r", "b"}, "units");
    const std::pair<const char*, const char*> expected[] = {{"m_prime", "1"}, {"theta", "rad"}, {"m", "m"}, {"r", "m"}, {"b", "1"}};
    for (const auto& [key, unit] : expected)
        if (text(units, key, "units") != unit)
            fail(std::format("units: '{}' is '{}', version 1 needs '{}'", key, units.at(key).get<std::string>(), unit));
    GeneralMesh mesh;
    if (root.contains("description")) mesh.description = text(root, "description", "the file");
    mesh.pitch = real(need(root, "pitch", "the file"), "pitch");
    const json& blocks = need(root, "blocks", "the file");
    if (!blocks.is_array()) fail("'blocks' must be an array");
    for (std::size_t k = 0; k < blocks.size(); ++k) {
        const std::string where = std::format("block {}", k);
        const json& b = blocks[k];
        only_keys(b, {"name", "ni", "nj", "m_prime", "theta", "m", "r", "b"}, where);
        MeshBlock B;
        if (b.contains("name")) B.name = text(b, "name", where);
        B.ni = integer(b, "ni", where);
        B.nj = integer(b, "nj", where);
        B.x = reals(b, "m_prime", where);
        B.y = reals(b, "theta", where);
        B.m = reals(b, "m", where);
        B.r = reals(b, "r", where);
        B.b = reals(b, "b", where);
        mesh.blocks.push_back(std::move(B));
    }
    const json& segments = need(root, "segments", "the file");
    if (!segments.is_array()) fail("'segments' must be an array");
    for (std::size_t k = 0; k < segments.size(); ++k) {
        const std::string where = std::format("segment {}", k);
        const json& s = segments[k];
        only_keys(s, {"block", "face", "range", "kind", "blade", "surface"}, where);
        FaceSegment g;
        g.block = integer(s, "block", where);
        g.face = face_of(text(s, "face", where), where);
        const json& range = need(s, "range", where);
        if (!(range.is_array() && range.size() == 2 && range[0].is_number_integer() && range[1].is_number_integer()))
            fail(where + ": 'range' must be [begin, end] (integers)");
        g.begin = range[0].get<int>();
        g.end = range[1].get<int>();
        const std::string kind = text(s, "kind", where);
        bool known = false;
        for (SegmentKind c : {SegmentKind::Blade, SegmentKind::Inlet, SegmentKind::Outlet, SegmentKind::Periodic,
                              SegmentKind::Interface})
            if (kind == segment_kind_name(c)) {
                g.kind = c;
                known = true;
            }
        if (!known) fail(std::format("{}: kind '{}' is not blade, inlet, outlet, periodic or interface", where, kind));
        if (g.kind == SegmentKind::Blade) {
            g.blade = integer(s, "blade", where);
            g.surface = surface_of(text(s, "surface", where), where);
        } else if (s.contains("blade") || s.contains("surface")) {
            fail(where + ": only a blade segment has 'blade' and 'surface'");
        }
        mesh.segments.push_back(g);
    }
    const json& connections = need(root, "connections", "the file");
    if (!connections.is_array()) fail("'connections' must be an array");
    for (std::size_t k = 0; k < connections.size(); ++k) {
        const std::string where = std::format("connection {}", k);
        const json& c = connections[k];
        only_keys(c, {"kind", "segments", "reversed", "image"}, where);
        const std::string kind = text(c, "kind", where);
        if (kind != "periodic" && kind != "interface") fail(std::format("{}: kind '{}' is not periodic or interface", where, kind));
        const json& pair = need(c, "segments", where);
        if (!(pair.is_array() && pair.size() == 2 && pair[0].is_number_integer() && pair[1].is_number_integer()))
            fail(where + ": 'segments' must be [first, second] (segment indices)");
        Connection cn;
        cn.first = pair[0].get<int>();
        cn.second = pair[1].get<int>();
        const json& rev = need(c, "reversed", where);
        if (!rev.is_boolean()) fail(where + ": 'reversed' must be true or false");
        cn.reversed = rev.get<bool>();
        if (kind == "periodic") {
            cn.image = integer(c, "image", where);
        } else if (c.contains("image")) {
            fail(where + ": an interface has no 'image'");
        }
        const int ns = static_cast<int>(mesh.segments.size());
        for (int s : {cn.first, cn.second})
            if (s < 0 || s >= ns) fail(std::format("{} names segment {}, the file has {}", where, s, ns));
        const SegmentKind want = kind == "periodic" ? SegmentKind::Periodic : SegmentKind::Interface;
        for (int s : {cn.first, cn.second})
            if (mesh.segments[static_cast<std::size_t>(s)].kind != want)
                fail(std::format("{} is a {} connection but segment {} is {}", where, kind, s,
                                 segment_kind_name(mesh.segments[static_cast<std::size_t>(s)].kind)));
        mesh.connections.push_back(cn);
    }
    const json& tags = need(root, "tags", "the file");
    if (!tags.is_array()) fail("'tags' must be an array");
    for (std::size_t k = 0; k < tags.size(); ++k) {
        const std::string where = std::format("tag {}", k);
        const json& t = tags[k];
        only_keys(t, {"edge", "blade", "surface", "block", "node"}, where);
        EdgeTag e;
        const std::string edge = text(t, "edge", where);
        if (edge == "leading") e.edge = EdgeKind::Leading;
        else if (edge == "trailing") e.edge = EdgeKind::Trailing;
        else fail(std::format("{}: edge '{}' is not leading or trailing", where, edge));
        e.blade = integer(t, "blade", where);
        e.surface = surface_of(text(t, "surface", where), where);
        e.block = integer(t, "block", where);
        const json& node = need(t, "node", where);
        if (!(node.is_array() && node.size() == 2 && node[0].is_number_integer() && node[1].is_number_integer()))
            fail(where + ": 'node' must be [i, j] (integers)");
        e.i = node[0].get<int>();
        e.j = node[1].get<int>();
        mesh.tags.push_back(e);
    }
    validate_mesh(mesh);
    return mesh;
}

void write_mesh_file(const GeneralMesh& mesh, const std::string& path)
{
    const std::string s = write_mesh_json(mesh);
    std::ofstream out(path, std::ios::binary);
    if (!out) fail(std::format("cannot open '{}' for writing", path));
    out << s;
    if (!out) fail(std::format("writing '{}' failed", path));
}

GeneralMesh read_mesh_file(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) fail(std::format("cannot open '{}'", path));
    std::ostringstream s;
    s << in.rdbuf();
    return read_mesh_json(s.str());
}

}  // namespace BladeGen::Turbolenza
