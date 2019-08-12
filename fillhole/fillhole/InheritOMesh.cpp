#include "InheritOMesh.h"
#include <cstdio>
#include <string>
#include <vector>
#include <queue>

#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>
//#include <OpenMesh/Tools/Subdivider/Uniform/LongestEdgeT.hh>// dead loop?
#include <OpenMesh/Tools/Subdivider/Uniform/Sqrt3T.hh>
//#include <OpenMesh/Tools/Subdivider/Uniform/Sqrt3InterpolatingSubdividerLabsikGreinerT.hh>//error
//#include <OpenMesh/Tools/Subdivider/Uniform/ModifiedButterFlyT.hh>////error
//#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>//quad poly
#include <OpenMesh/Tools/Smoother/SmootherT.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

using namespace std;

#define gydebug

#ifdef gydebug
static int debugCnt = 0;
#endif // gydebug



bool compareVertexAngleAsc(const VH & small, const VH & big)
{
	//return mesh.property(vertexAngle, big) > mesh.property(vertexAngle, small);//big .angle > small.angle
	return big.edgeAngle > small.edgeAngle;
}

int MeshHandle::ReadMeshFile(string file) {
	// 请求顶点法线 vertex normals
	mesh.request_vertex_normals();
	//如果不存在顶点法线，则报错 
	if (!mesh.has_vertex_normals())
	{
		cout << "错误：标准定点属性 “法线”不存在" << endl;
		return -1;
	}

	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, file, opt))
	{
		cout << "无法读取文件:" << file << endl;
		return -1;
	}

	//如果不存在顶点法线，则计算出
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// 通过面法线计算顶点法线
		mesh.request_face_normals();
		// mesh计算出顶点法线
		mesh.update_normals();
		// 释放面法线
		//mesh.release_face_normals();
	}

	std::cout << "ReadMeshFile done" << file << std::endl;
	return 0;
}

int MeshHandle::Init(string infile) {

	ReadMeshFile(infile);
	NeighborRingNums = 3;

	std::cout << "Init done" << std::endl;
	return 0;
}

int MeshHandle::WriteMeshFile(std::string file)
{
	OpenMesh::IO::write_mesh(mesh, file);
	std::cout << "WriteMeshFile done" << file << std::endl;
	return 0;
}

int MeshHandle::Save_file_as(std::string infile, std::string outfile)
{
	return ReadMeshFile(infile) || WriteMeshFile(outfile);
}

int MeshHandle::FindHoles()
{
	//OpenMesh::VPropHandleT<int> vertexType;
	//mesh.add_property(vertexType);
	holes.clear();
	for (auto vertexIt = mesh.vertices_begin(); vertexIt != mesh.vertices_end(); vertexIt++) 
	{
		//find a boundary vertex, means a hole start here 
		if ( mesh.is_boundary(*vertexIt) && mesh.data(*vertexIt).type < 1 )//mesh.property(vertexType, *vertexIt) < 1)
		{
			//mesh.property(vertexType, *vertexIt) = 1;
			mesh.data(*vertexIt).type = 1;
			SingleHole aHole;
			aHole.AvergEdgeLen = 0;
			aHole.StartVertex = *vertexIt;
			aHole.VertBoundaryHs.push_back(*vertexIt);
			aHole.CenterPoint = mesh.point(*vertexIt);
			OpenMesh::VertexHandle nextVertex= *vertexIt;
			//along boundary push all boundary vertex
			do{
				//find a boundary neighbor vertex
				for (auto neighborIt = mesh.cvv_begin(nextVertex); neighborIt != mesh.cvv_end(nextVertex); neighborIt++)
				{
					if (*neighborIt == aHole.StartVertex)
					{
						nextVertex = *neighborIt; //break while
						break;
					}
					if (mesh.is_boundary(*neighborIt) && mesh.data(*neighborIt).type < 1)// mesh.property(vertexType, *neighborIt) < 1) //  wrong: vertex with more than one hole
					{
						//mesh.property(vertexType, *neighborIt) = 1;
						mesh.data(*neighborIt).type = 1;
						aHole.VertBoundaryHs.push_back(*neighborIt);
						aHole.CenterPoint += mesh.point(*neighborIt);
						aHole.AvergEdgeLen += ( mesh.point(nextVertex) - mesh.point(*neighborIt) ).length();
						nextVertex = *neighborIt;
						break;// assumpt only one  boundary vertex, so got find next
					}
					//else if (0 == mesh.data(*neighborIt).type) //neighbor ring 2
					//{
					//	mesh.data(*neighborIt).type = 2; // push to a queue;
					//	//debug
					//	mesh.set_point(*neighborIt, mesh.point(*neighborIt) + mesh.normal(*neighborIt) );
					//}
				}
				//
			}while (nextVertex != aHole.StartVertex && nextVertex.is_valid());

			if (aHole.VertBoundaryHs.size() < 3) {
				continue;//false hole
			}
			aHole.AvergEdgeLen /= (aHole.VertBoundaryHs.size() - 1);
			aHole.CenterPoint /= aHole.VertBoundaryHs.size();

			holes.push_back(aHole);
		}
	}
	return 0;
}

/*
pre processing for example: split long edge and  collapse  short edge ,(0.5-2) or (0.8-1.33 )

有尖刺穿透面片 FillHoleMinAngle 最后一个洞可能是尖刺（倒数第二个洞，四边形是凹），中间会形成很多尖刺

点的两边夹角应该加pi，当已经有一条边链接前后两个点，还有更糟糕的，总是要处理夹角大雨180

要考虑两个边的面法向尽可能平行，例如补一个锐角棱的洞，应该想办法在两个平面上假面，最后加棱
*/
int MeshHandle::FillHoles()
{
	cout << "FillHoles..." << endl;
	for (vector<SingleHole>::iterator h = holes.begin(); h != holes.end(); h++)
	{
		//judge
		//cout << "FillHoles...1" << endl;
		// pre processing for example: split long edge and  collapse  short edge ,(0.5-2) or (0.8-1.33 ) 
		//FillHoleCenterPoint(*h);// or sort area
		FillHoleMinAngle(*h);
		//SubdiveAHole(*h);
		//SmoothAHole(*h);//useless
	}
	//mesh.update_normals();
	////subdivision
	//OpenMesh::Subdivider::Uniform::LoopT<OMTriMesh> loopSubdivider;//FillHoleCenterPoint=> harsh hole
	//loopSubdivider.attach(mesh);
	//loopSubdivider(2);
	//loopSubdivider.detach();

	//OpenMesh::Subdivider::Uniform::Sqrt3T<OMTriMesh> subdivider;//FillHoleCenterPoint=> harsh hole
	//subdivider.attach(mesh);
	//subdivider(2);
	//subdivider.detach();


	//// smoother
	//OpenMesh::Smoother::JacobiLaplaceSmootherT<OMTriMesh> smoother(mesh);// FillHoleMinAngle=>only leftharsh mesh
	//smoother.initialize(OpenMesh::Smoother::SmootherT<OMTriMesh>::Component::Tangential_and_Normal,   //Smooth direction                    
	//	OpenMesh::Smoother::SmootherT<OMTriMesh>::Continuity::C2);                      //Continuity 
	//smoother.smooth(3);

	cout << "FillHoles...end" << endl;
	return 0;
}

int MeshHandle::FillHoleCenterPoint(SingleHole & hole)
{
	//it is sure hole.VertHssize >= 3
	OpenMesh::VertexHandle CenterV = mesh.add_vertex(hole.CenterPoint);
	mesh.data(CenterV).type = -1;//new point
	auto vIt = hole.VertBoundaryHs.begin();
	auto vItNext = vIt+1;
	for (; vItNext != hole.VertBoundaryHs.end(); vIt++, vItNext++)
	{
		//cout << "debug... a " << vIt->idx() << " b " << vItNext->idx() << " c " << CenterV.idx() << endl;
		hole.NewtFaceHs.push_back( mesh.add_face(*vIt, *vItNext, CenterV) );//must be clock wise
	}
	hole.NewtFaceHs.push_back( mesh.add_face(*vIt, hole.StartVertex , CenterV) );
	return !( hole.VertBoundaryHs.size() == hole.NewtFaceHs.size() );
}

int MeshHandle::FillHoleMinAngle(SingleHole & hole)
{
	if (3 > hole.VertBoundaryHs.size()) {
		cerr << "false hole" << endl;
		return -1;
	}

	vector<VH> tempVs;
	calculateVerteciesAngle(tempVs, hole);
	std::vector<VH>::iterator preVH,minVH,nextVH;
	while (tempVs.size() > 3 ) { // debug ball has 9 vertex
		minVH = std::min_element(std::begin(tempVs), std::end(tempVs), [](const VH & a, const VH & b)->bool {
			return a.edgeAngle < b.edgeAngle;
		});//compareVertexAngleAsc

		try {
			//todo 
			GetPrevNextIter(tempVs, preVH, minVH, nextVH);
			OpenMesh::FaceHandle newF = mesh.add_face(preVH->vh, minVH->vh, nextVH->vh);
			if (newF.is_valid())
			{
				hole.NewtFaceHs.push_back(newF);
			}
			//else debug
			//tempVs remove this min VH
			deleteMinVertex(tempVs, minVH);
		}
		catch (std::exception& x) {
			std::cerr << x.what() << std::endl;
			WriteMeshFile("G:\\gyGit\\data\\hole\\out_temp_err.obj");
		}
#ifdef gydebug
		//WriteMeshFile("G:\\gyGit\\data\\hole\\out_temp_d"+ to_string(debugCnt++) +".obj");
#endif // gydebug
	}
	//todo watch out ! it colude be concave
	if (3== tempVs.size()) {
		try{
			// are angles right?
			//adjustOneConcaveVertex(tempVs[0].vh, tempVs[1], tempVs[2].vh);
			//adjustOneConcaveVertex(tempVs[1].vh, tempVs[2], tempVs[0].vh);
			//adjustOneConcaveVertex(tempVs[2].vh, tempVs[1], tempVs[0].vh);

			OpenMesh::FaceHandle newF = mesh.add_face(tempVs[0].vh, tempVs[1].vh, tempVs[2].vh);
			hole.NewtFaceHs.push_back(newF);
		}
		catch (std::exception& x) {
			std::cerr << x.what() << std::endl;
		}
	}
	else {
		cerr << "fill error,still hole" << endl;
		return -1;
	}
	return 0;
}


int MeshHandle::SubdiveAHole(SingleHole & hole)// refine a bad triangle
{
	bool isNeedSubdive = true;
	vector<OpenMesh::FaceHandle> facesH;
	facesH.clear();
	flipFaces(hole.NewtFaceHs);
	mesh.request_face_status();
	
	while (isNeedSubdive) {
		isNeedSubdive = false;
		for (vector<OpenMesh::FaceHandle>::iterator fIt = hole.NewtFaceHs.begin(); fIt != hole.NewtFaceHs.end(); fIt++) {
			//each face if longest > shortest * 3 need subdive or flip
			OpenMesh::Vec3f barycenter;
			OpenMesh::VertexHandle pointA, pointB, pointC;

			auto vIt = mesh.fv_ccwbegin(*fIt);
			pointA = *vIt;
			pointB = *(++vIt);
			pointC = *(++vIt);
			barycenter = mesh.point(pointA) + mesh.point(pointB) + mesh.point(pointC);
			if (++vIt != mesh.fv_ccwend(*fIt)) {
				continue;
			}
			//int vertexNums = 3;// it must be three //todo judge vIt size
			barycenter /= 3;
			if ((mesh.point(pointB) - mesh.point(pointA)).length() > hole.AvergEdgeLen
				|| (mesh.point(pointC) - mesh.point(pointB)).length() > hole.AvergEdgeLen
				|| (mesh.point(pointA) - mesh.point(pointC)).length() > hole.AvergEdgeLen)
			{

				OpenMesh::VertexHandle newV = mesh.add_vertex(barycenter);
				mesh.data(newV).type = -1;//new point
				if (!mesh.status(*fIt).deleted()) {
					try {
						//mesh.split();
						isNeedSubdive = true;
						mesh.delete_face(*fIt, false);
						OpenMesh::FaceHandle fa = mesh.add_face(newV, pointA, pointB);
						OpenMesh::FaceHandle fb = mesh.add_face(newV, pointB, pointC);
						OpenMesh::FaceHandle fc = mesh.add_face(newV, pointC, pointA);

						//add new face to arrary lisr
						//hole.NewtFaceHs.erase(fIt);
						if (fa.is_valid() && fb.is_valid() && fc.is_valid()) {
							facesH.push_back(fa);
							facesH.push_back(fb);
							facesH.push_back(fc);
						}

					}
					catch (std::exception& x) {
						std::cerr << x.what() << std::endl;
					}
				}
			}
			else {
				//need not subdive
				facesH.push_back(*fIt);
			}

		}
		flipFaces(facesH);
	}
	hole.NewtFaceHs.clear();
	hole.NewtFaceHs = facesH;
	//mesh.garbage_collection();//NOTICE this operation order
	return 0;
}

int MeshHandle::SmoothAHole(SingleHole & hole)
{
	OpenMesh::Vec3f barycenter(0, 0, 0);
	int nums = 0;
	for (vector<OpenMesh::FaceHandle>::iterator fIt = hole.NewtFaceHs.begin(); fIt != hole.NewtFaceHs.end(); fIt++) {
		for (auto vIt = mesh.fv_cwbegin(*fIt); vIt != mesh.fv_cwend(*fIt); vIt++)
		{
			if (mesh.data(*vIt).type <0 ) {
				barycenter = OpenMesh::Vec3f(0, 0, 0);
				nums = 0;
				for (auto vvIt = mesh.vv_cwbegin(*vIt); vvIt != mesh.vv_cwend(*vIt); vvIt++) {
					barycenter += mesh.point(*vIt);
					nums++;
				}
				barycenter /= nums;
				if (nums> 100) {
					cout<<"god!"<<nums<<endl;
				}
				mesh.set_point(*vIt,barycenter);
			}
			
		}

	}
	return 0;
}

int MeshHandle::deleteMinVertex(std::vector<VH>& tempVs, std::vector<VH>::iterator minVH)
{
	std::vector<VH>::iterator leftVH, preVH, nextVH, rightVH;

	preVH = GetPrevIter(tempVs,minVH);
	tempVs.erase(minVH);

	//after delete minvh point t o next so
	//re find nextVH is next 
	//and preVH is need update for get left of preVH
	nextVH = GetNextIter(tempVs, preVH);
	preVH = GetPrevIter(tempVs, nextVH);
	leftVH = GetPrevIter(tempVs, preVH);
	rightVH = GetNextIter(tempVs, nextVH);

	//preVH->orgArrIndexNext = nextVH->orgArrIndex;
	//nextVH->orgArrIndexPre = preVH->orgArrIndex;
#ifdef gydebug
	//if (preVH->orgArrIndex == 22) {
	//	cout << "l" << preVH->toString() << endl;
	//	cout << "r" << nextVH->toString() << endl << endl;
	//}
#endif // gydebug
	calculateOneVertexAngle(leftVH->vh, preVH, nextVH->vh);
	calculateOneVertexAngle(preVH->vh, nextVH, rightVH->vh);
	return 0;
}

int MeshHandle::calculateOneVertexAngle(OpenMesh::VertexHandle pre, VH & cur, OpenMesh::VertexHandle next)
{
	//param must be valid
	OpenMesh::Vec3f lineA, lineB;
	lineA = mesh.point(pre) - mesh.point(cur.vh);
	lineB = mesh.point(next) - mesh.point(cur.vh);

	OpenMesh::Vec3f faceNewNormal, facePreNormal, facenextNormal;
	facePreNormal = mesh.calc_face_normal(mesh.face_handle(mesh.find_halfedge(cur.vh, pre)));//from vstart to vend
	facenextNormal = mesh.calc_face_normal(mesh.face_handle(mesh.find_halfedge(next, cur.vh)));

	cur.edgeAngle = acos(OpenMesh::dot(lineA.normalize(), lineB.normalize()));
	cur.edgeFacesNormalAngle = acos(OpenMesh::dot(facePreNormal, facenextNormal));

	if (mesh.find_halfedge(pre, next).is_valid() )  { //  edge exit or new face normal donnot parallel
		//|| M_PI / 2 > acos( OpenMesh::dot( faceNewNormal, ( facePreNormal + facenextNormal ) / 2 )
		cur.edgeAngle += M_PI; // obtuse angle
	}
	else {
		faceNewNormal = OpenMesh::cross(lineA, lineB);//notice the dirct
		//faceNewNormal = mesh.calc_face_normal(mesh.point(pre), mesh.point(cur), mesh.point(next));//OpenMesh::cross();
		if ( M_PI / 2 > acos(OpenMesh::dot(faceNewNormal, (facePreNormal + facenextNormal) / 2) ) ) {
			cur.edgeAngle += M_PI; // obtuse angle
		}
	}
	return 0;
}

int MeshHandle::calculateOneVertexAngle(OpenMesh::VertexHandle pre, std::vector<VH>::iterator & cur, OpenMesh::VertexHandle next)
{
	OpenMesh::Vec3f lineA, lineB;
	lineA = mesh.point(pre) - mesh.point(cur->vh);
	lineB = mesh.point(next) - mesh.point(cur->vh);

	OpenMesh::Vec3f faceNewNormal, facePreNormal, facenextNormal;
	facePreNormal = mesh.calc_face_normal(mesh.face_handle(mesh.find_halfedge(cur->vh, pre)));//from vstart to vend
	facenextNormal = mesh.calc_face_normal(mesh.face_handle(mesh.find_halfedge(next, cur->vh)));

	cur->edgeAngle = acos(OpenMesh::dot(lineA.normalize(), lineB.normalize()));
	cur->edgeFacesNormalAngle = acos(OpenMesh::dot(facePreNormal, facenextNormal));

	if (mesh.find_halfedge(pre, next).is_valid()) { //  edge exit or new face normal donnot parallel
		//|| M_PI / 2 > acos( OpenMesh::dot( faceNewNormal, ( facePreNormal + facenextNormal ) / 2 )
		cur->edgeAngle += M_PI; // obtuse angle
	}
	else {
		faceNewNormal = OpenMesh::cross(lineA, lineB);//notice the dirct
		//faceNewNormal = mesh.calc_face_normal(mesh.point(pre), mesh.point(cur), mesh.point(next));//OpenMesh::cross();
		if (M_PI / 2 > acos(OpenMesh::dot(faceNewNormal, (facePreNormal + facenextNormal) / 2))) {
			cur->edgeAngle += M_PI; // obtuse angle
		}
	}
	return 0;
}

int MeshHandle::calculateVerteciesAngle(std::vector<VH> & tempVs, SingleHole & hole)
{
	tempVs.clear();
	//int i = 1;
	auto vItPrevious = hole.VertBoundaryHs.begin();
	auto vIt = vItPrevious + 1;
	auto vItNext = vIt + 1;

	//fisrt point is vItPrevious
	VH vfirst;
	vfirst.vh = (*vItPrevious);
	//vfirst.orgArrIndex = 0;
	//vfirst.orgArrIndexPre = hole.VertBoundaryHs.size() - 1;
	//	vfirst.orgArrIndexNext = 1;
	tempVs.push_back(vfirst);
	//for 1 : n-1
	for (; vItNext != hole.VertBoundaryHs.end(); vItPrevious++, vIt++, vItNext++)
	{
		VH v;
		v.vh = (*vIt);
		//v.orgArrIndex = i;
		//v.orgArrIndexPre = i - 1;
		//v.orgArrIndexNext = i + 1;
		//i++
		calculateOneVertexAngle(*vItPrevious, v, *vItNext); 
		tempVs.push_back(v);
	}
	//assert i == hole.VertBoundaryHs.size()-1
	//last point is vIt
	VH vlast;
	vlast.vh = (*vIt);
	//vlast.orgArrIndex = i; // size()-1
	//vlast.orgArrIndexPre = i - 1;
	//vlast.orgArrIndexNext = 0;
	calculateOneVertexAngle(*vItPrevious, vlast, tempVs[0].vh); //Vn-1 <- Vn -> V0
	tempVs.push_back(vlast);
	//fisrt point	
	calculateOneVertexAngle(*vIt, tempVs[0], tempVs[1].vh); //Vn <- V0 -> V1

	return 0;
}

int MeshHandle::flipFaces(std::vector<OpenMesh::FaceHandle>& faces)
{
	//next to do sort long edge and flip
	for (vector<OpenMesh::FaceHandle>::iterator fIt = faces.begin(); fIt != faces.end(); fIt++) {

		//flip
		OpenMesh::EdgeHandle heh[4];
		int i = 0;
		for (auto hIt = mesh.fh_cwbegin(*fIt); hIt != mesh.fh_cwend(*fIt); hIt++) {
			heh[i++] = mesh.edge_handle(*hIt);
			//if (mesh.is_flip_ok(mesh.edge_handle(*hIt))) {
			//	//mesh.flip( mesh.edge_handle(*hIt) );
			//}
		}

		heh[3] = mesh.calc_edge_length(heh[1]) > mesh.calc_edge_length(heh[0]) ? heh[1] : heh[0];
		heh[3] = mesh.calc_edge_length(heh[3]) > mesh.calc_edge_length(heh[2]) ? heh[3] : heh[2];
		
		if (mesh.is_flip_ok(heh[3])) {
			mesh.flip(heh[3]);
		}
	}

	return 0;
}

int MeshHandle::adjustOneConcaveVertex(OpenMesh::VertexHandle pre, VH & cur, OpenMesh::VertexHandle next)
{
	//adjust
	if (M_PI > cur.edgeAngle) { // 180 145 120 90 60 
		return 0;//need not
	}
	OpenMesh::Vec3f temp(0,0,0);

	//todo
	// maybe there is not other points who is not boundary,so remove this is best option? ==> if (0 == nums)
	// maybe other neighbor points are aslo in the hole side of line between  pre next
	int nums = 0;
	for (auto vvit = mesh.vv_cwbegin(cur.vh); vvit != mesh.vv_cwend(cur.vh);vvit++, nums++) {
		//if not boundary
		temp += mesh.point(*vvit);
	}
	temp -= mesh.point(pre);
	temp -= mesh.point(next);
	nums -= 2;

	if (0 == nums) {
		return -1;
	}
	temp += mesh.point(cur.vh);
	temp /= (nums + 1);

	mesh.set_point(cur.vh, temp);
	calculateOneVertexAngle(pre, cur, next);

	return 0;
}



///////////////////////////////////////////////////////////////////////
