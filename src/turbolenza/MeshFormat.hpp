/// @file turbolenza/MeshFormat.hpp
/// @brief The "turbolenza.mesh" file format, version 1: a GeneralMesh as JSON (the contract with
///        ENIGMAgrid's writer; specification docs/turbolenza/mesh_format.md).  Not a port of v1.
///
/// The writer prints every number with 17 significant digits ({:.17g}: a double reads back bit for bit)
/// and states the units (m' dimensionless, theta rad, m and r in metres, b relative).  The reader is
/// strict: the format name, version 1, the units, every required key with its type, no unknown key;
/// then validate_mesh.  Anything else is refused with a MeshFormatError (schema) or a MeshError (an
/// inconsistent mesh: unpaired periodic segment, interface node-count mismatch, negative Jacobian,
/// missing blade tags, ...).  Dirichlet segments (verification only) are not part of the format.
#pragma once

#include "turbolenza/GeneralMesh.hpp"

#include <string>

namespace BladeGen::Turbolenza {

inline constexpr const char* kMeshFormatName = "turbolenza.mesh";
inline constexpr int kMeshFormatVersion = 1;

class MeshFormatError : public MeshError
{
public:
    using MeshError::MeshError;
};

/// The JSON text of a valid mesh (validate_mesh first; MeshFormatError for a Dirichlet segment).
[[nodiscard]] std::string write_mesh_json(const GeneralMesh& mesh);
/// The mesh of a JSON text (see the file comment).
[[nodiscard]] GeneralMesh read_mesh_json(const std::string& text);

/// The same through a file (MeshFormatError when it cannot be written or read).
void write_mesh_file(const GeneralMesh& mesh, const std::string& path);
[[nodiscard]] GeneralMesh read_mesh_file(const std::string& path);

}  // namespace BladeGen::Turbolenza
