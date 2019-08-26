#pragma once


#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriConnectivity.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh/Core/Mesh/ArrayKernel.hh>



struct MyTraits :public OpenMesh::DefaultTraits {
public:
	VertexTraits{
		//Point cog; // center of gravity重心
		int type=0;// 1:boundray -1:added_vertex
	};
	FaceTraits{
		//Point cog; // center of gravity重心
	};
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> OMTriMesh;


class VH {
public:
	OpenMesh::VertexHandle vh;
	//int orgArrIndexPre;
	//int orgArrIndex;
	//int orgArrIndexNext;
	double edgeAngle;// 0-360; must < 180
	//todo use it 
	double edgeFacesNormalAngle;//0-180 ; should <30 or less

	std::string toString() {
		std::string ret;
		ret = " vh "+ std::to_string(vh.idx())
			+ "\n edgeAngle " + std::to_string(edgeAngle * 180 / M_PI) + " f " + std::to_string(edgeFacesNormalAngle*180/M_PI)
			//+ "\n p " + std::to_string(orgArrIndexPre) + " o " + std::to_string(orgArrIndex) + " n " + std::to_string(orgArrIndexNext)
			;
		return ret;
	}
};

//VertexHandle :single hole of mesh
class SingleHole {
public:
	//int NeighborRingNums;
	double AvergEdgeLen;
	OpenMesh::Vec3f CenterPoint;
	OpenMesh::VertexHandle StartVertex;
	std::vector<OpenMesh::VertexHandle> VertBoundaryHs;//boundary of origin mesh

	std::vector<OpenMesh::FaceHandle> NewtFaceHs;
	//std::map<int, int> VertexType;
	//判断类型，是否是一个凹的多边形，类似U型，中心在洞外，不可以中心补洞
};


//mesh with hole depose method
class MeshHandle {
public:	
	OMTriMesh mesh;
	
	std::vector<SingleHole> holes;
	int NeighborRingNums;

	int Init(std::string infile);

	int ReadMeshFile(std::string file);
	int WriteMeshFile(std::string file);

	int Save_file_as(std::string infile, std::string outfile);
	
	int FindHoles();

	int FillHoles();
	int FillHoleCenterPoint(SingleHole & hole);
	int FillHoleMinAngle(SingleHole & hole);// WithSplitLongEde
	// FillHoleDihedralAngle(SingleHole & hole);
	// min mesh.calc_dihedral_angle
protected:
private:
	int SubdiveAHole(SingleHole & hole);
	int SmoothAHole(SingleHole & hole);

	//OpenMesh::VPropHandleT<double> vertexAngle;
	int deleteMinVertex(std::vector<VH> & tempVs, std::vector<VH>::iterator minVH);
	int calculateOneVertexAngle(OpenMesh::VertexHandle pre, VH & cur, OpenMesh::VertexHandle next);//isObtuseAngle();
	int calculateOneVertexAngle(OpenMesh::VertexHandle pre, std::vector<VH>::iterator & cur, OpenMesh::VertexHandle next);//isObtuseAngle();
	int calculateVerteciesAngle(std::vector<VH> & tempVs, SingleHole & hole);
	int flipFaces(std::vector<OpenMesh::FaceHandle> & faces);
	//bool compareVertexAngleAsc(const VH & big, const VH & small);

	int adjustOneConcaveVertex(OpenMesh::VertexHandle pre, VH & cur, OpenMesh::VertexHandle next);
};
///////////////////////////////////////////////////////////////////////
template <class TVH>
int GetPrevNextIter(const std::vector<TVH> & tempVs, typename std::vector<TVH>::const_iterator & prevVH, typename std::vector<TVH>::iterator cur, typename std::vector<TVH>::const_iterator & nextVH)
{
	//cur = std::find(std::begin(tempVs), std::end(tempVs), *this);
	if (cur != tempVs.begin()) {
		prevVH = std::prev(cur);
		nextVH = std::next(cur);
		if (cur == tempVs.end() || nextVH == tempVs.end()) {
			nextVH = tempVs.begin();
		}
	}
	else {//minVH is begin the first one
		prevVH = std::prev(tempVs.end());
		nextVH = std::next(cur);// nextVH is second one,so must not end 
	}
	return 0;
}

template <class TVH>
typename std::vector<TVH>::iterator   GetPrevIter(std::vector<TVH> & tempVs, typename std::vector<TVH>::iterator cur)
{
	if (cur == tempVs.begin()) {
		return  std::prev(tempVs.end());
	}
	return std::prev(cur);
}

template <class TVH>
typename std::vector<TVH>::iterator   GetNextIter(std::vector<TVH> & tempVs,  typename std::vector<TVH>::iterator cur)
{
	if (cur == tempVs.end() || std::next(cur) == tempVs.end()) {
		return tempVs.begin();
	}
	return std::next(cur);
}
//////////////////////////////////////////////////////////////////////////
