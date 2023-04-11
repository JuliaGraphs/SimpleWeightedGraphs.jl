"""
    _get_nz_index!(mat::SparseMatrixCSC, i, j)

Return the index in nzval of `mat[i, j]`. We assume bounds are already checked

See https://github.com/JuliaSparse/SparseArrays.jl/blob/fa547689947fadd6c2f3d09ddfcb5f26536f18c8/src/sparsematrix.jl#L2492 for implementation
"""
function _get_nz_index!(mat::SparseMatrixCSC, i::Integer, j::Integer)
    # r1 and r2 are start and end of the column
    @inbounds r1 = Int(mat.colptr[j])
    @inbounds r2 = Int(mat.colptr[j + 1] - 1)
    (r1 > r2) && return 0 # column is empty so we have a non structural zero
    # search if i correspond to a stored value
    @inbounds indx = searchsortedfirst(mat.rowval, i, r1, r2, Base.Forward)
    @inbounds ((indx > r2) || (mat.rowval[indx] != i)) && return 0
    return indx
end
