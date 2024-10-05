#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <iostream>
#include <string>
#include <list>
#include <fstream>
#include <cstring>
#include <math.h>
#include <vector>

int min_treshold = 0.01f;

typedef struct ponto {
    float x;
    float y;
    float z;
} *PONTO;

typedef struct objeto {
    std::vector<PONTO> vertices;
    std::vector<unsigned int> indices;
    std::vector<PONTO> tex_vert;
    std::vector<PONTO> normals;
} *OBJETO;


float bezierM[4][4] = {{ -1.0f, 3.0f, -3.0f, 1.0f},
					   { 3.0f, -6.0f, 3.0f, 0.0f},
					   { -3.0f, 3.0f, 0.0f, 0.0f},
					   { 1.0f, 0.0f, 0.0f, 0.0f}};


void multMatrixVector(float *m, float *v, float *res) {

	for (int j = 0; j < 4; ++j) {
		res[j] = 0;
		for (int k = 0; k < 4; ++k) {
			res[j] += v[k] * m[j * 4 + k];
		}
	}
}

void cross(float *a, float *b, float *res) {

	res[0] = a[1]*b[2] - a[2]*b[1];
	res[1] = a[2]*b[0] - a[0]*b[2];
	res[2] = a[0]*b[1] - a[1]*b[0];
}

void normalize(float *a) {

	float l = sqrt(a[0]*a[0] + a[1] * a[1] + a[2] * a[2]);
	a[0] = a[0]/l;
	a[1] = a[1]/l;
	a[2] = a[2]/l;
}

void normalize_normals(OBJETO obj){

    for (PONTO p : obj->normals){
        float c[3] = {p->x, p->y, p->z};
        normalize(c);
        p->x = c[0];
        p->y = c[1];
        p->z = c[2];
    }
}

bool compare_pontos(PONTO p1, PONTO p2){
    if(abs(p1->x - p2->x) < min_treshold &&
       abs(p1->y - p2->y) < min_treshold &&
       abs(p1->z - p2->z) < min_treshold)
            return true;
    return false;
}

int ponto_in_vector(PONTO p, std::vector<PONTO> vec){
    for(int i = 0; i < vec.size(); i++)
        if(compare_pontos(p, vec[i]))
            return i;
    return -1;
}

void index_vertices(OBJETO obj){

    normalize_normals(obj);

    std::vector<PONTO> new_vert;
    std::vector<PONTO> new_tex;
    std::vector<PONTO> new_norm;

    for (int i = 0; i < obj->vertices.size(); i++){
        PONTO p = obj->vertices[i];
        PONTO tp = obj->tex_vert[i];
        PONTO np = obj->normals[i];
        int idx = ponto_in_vector(p, new_vert);
        if (idx == -1 ||
           (idx != -1 &&
             (!compare_pontos(tp, new_tex[idx]) ||
             !compare_pontos(np, new_norm[idx])))){
            idx = new_vert.size();
            new_vert.push_back(p);
            new_tex.push_back(tp);
            new_norm.push_back(np);
        }
        obj->indices.push_back(idx);
    }
    obj->vertices = new_vert;
    obj->tex_vert = new_tex;
    obj->normals = new_norm;
}


// Função que escreve para ficheiro
void makeFile(std::string fich, OBJETO obj) {
    std::ofstream file(fich);
    file << obj->vertices.size() << std::endl;
    file << obj->indices.size() << std::endl;
    file << obj->tex_vert.size() << std::endl;
    file << obj->normals.size() << std::endl;

    // Write vertices
    for (PONTO it : obj->vertices){
        file << it->x << std::endl;
        file << it->y << std::endl;
        file << it->z << std::endl;
    }

    // Write indices
    for (unsigned int it : obj->indices){
        file << it << std::endl;
    }

    // Write texture vertices
    for (PONTO it : obj->tex_vert){
        file << it->x << std::endl;
        file << it->y << std::endl;
    }

    // Write normals
    for (PONTO it : obj->normals){
        file << it->x << std::endl;
        file << it->y << std::endl;
        file << it->z << std::endl;
    }

    file.close();
}

// Função que adiciona pontos a uma std::list
void addPonto (float x, float y, float z, OBJETO obj) {
    PONTO p = new ponto();
    p->x = x;
    p->y = z;//y;
    p->z = y;//z;
    obj->vertices.push_back(p);
}

// Função que adiciona pontos a uma std::list
void addTexPonto (float u, float v, OBJETO obj) {
    PONTO p = new ponto();
    p->x = u;
    p->y = v;
    p->z = 0;
    obj->tex_vert.push_back(p);
}

// Função que adiciona pontos a uma std::list
void addNormPonto (float x, float y, float z, OBJETO obj) {
    PONTO p = new ponto();
    p->x = x;
    p->y = z;//y;
    p->z = y;//z;
    obj->normals.push_back(p);
}

// Função do plano
void planeFunc(float tamanho, int div, std::string fich) {
    
    OBJETO obj = new objeto();
    float step = (float) tamanho/div;
    float u = 0.0;
    float v = 1.0;
    float u_step = 1.0 / float(div);
    float v_step = 1.0 / float(div);

    for (float x = -tamanho/2.0f; x < tamanho/2.0f; x+=step ) {
        v = 1.0;
        for(float y = -tamanho/2.0f; y < tamanho/2.0f; y+=step) {
            //Triangulo 1
            addPonto(x,y,0,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x+step,y+step,0,obj);
            addTexPonto(u + u_step, v - v_step, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x+step,y,0,obj);
            addTexPonto(u + u_step, v, obj);
            addNormPonto(0, 0, 1, obj);

            //Triangulo 2
            addPonto(x,y,0,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x,y+step,0,obj);
            addTexPonto(u, v - v_step, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x+step,y+step,0,obj);
            addTexPonto(u + u_step, v - v_step, obj);
            addNormPonto(0, 0, 1, obj);

            v -= v_step;
        }
        u += u_step;
    }

    index_vertices(obj);
    makeFile(fich, obj);
}


// Função do box
void boxFunc(float tamanho, int div, std::string fich) {
    
    OBJETO obj = new objeto();
    float count = (float) tamanho/div;
    float u = 0.0;
    float v = 1.0;
    float u_step = 1.0 / float(div);
    float v_step = 1.0 / float(div);

    // Base inferior
    for (float x = -tamanho/2.0f; x < tamanho/2.0f; x+=count ) {
        v = 0.0;
        for(float y = -tamanho/2.0f; y < tamanho/2.0f; y+=count) {
            //Triangulo 1
            addPonto(x,y,-tamanho/2.0f,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 0, -1, obj);
            addPonto(x+count,y,-tamanho/2.0f,obj);
            addTexPonto(u + u_step, v, obj);
            addNormPonto(0, 0, -1, obj);
            addPonto(x+count,y+count,-tamanho/2.0f,obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(0, 0, -1, obj);

            //Triangulo 2
            addPonto(x,y,-tamanho/2.0f,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 0, -1, obj);
            addPonto(x+count,y+count,-tamanho/2.0f,obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(0, 0, -1, obj);
            addPonto(x,y+count,-tamanho/2.0f,obj);
            addTexPonto(u, v + v_step, obj);
            addNormPonto(0, 0, -1, obj);
            
            v += v_step;
        }
        u += u_step;
    }

    u = 0.0;
    v = 1.0;

    // Base superior
    for (float x = -tamanho/2.0f; x < tamanho/2.0f; x+=count ) {
        v = 1.0;
        for(float y = -tamanho/2.0f; y < tamanho/2.0f; y+=count) {
            //Triangulo 1
            addPonto(x,y,tamanho/2.0f,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x+count,y+count,tamanho/2.0f,obj);
            addTexPonto(u + u_step, v - v_step, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x+count,y,tamanho/2.0f,obj);
            addTexPonto(u + u_step, v, obj);
            addNormPonto(0, 0, 1, obj);

            //Triangulo 2
            addPonto(x,y,tamanho/2.0f,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x,y+count,tamanho/2.0f,obj);
            addTexPonto(u, v - v_step, obj);
            addNormPonto(0, 0, 1, obj);
            addPonto(x+count,y+count,tamanho/2.0f,obj);
            addTexPonto(u + u_step, v - v_step, obj);
            addNormPonto(0, 0, 1, obj);

            v -= v_step;
        }
        u += u_step;
    }

    u = 1.0;
    v = 0.0;

    // Base de trás
    for (float x = -tamanho/2.0f; x < tamanho/2.0f; x+=count ) {
        v = 0.0;
        for(float z = -tamanho/2.0f; z < tamanho/2.0f; z+=count) {
            //Triangulo 1
            addPonto(x,-tamanho/2.0f,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, -1, 0, obj);
            addPonto(x+count,-tamanho/2.0f,z+count,obj);
            addTexPonto(u - u_step, v + v_step, obj);
            addNormPonto(0, -1, 0, obj);
            addPonto(x+count,-tamanho/2.0f,z,obj);
            addTexPonto(u - u_step, v, obj);
            addNormPonto(0, -1, 0, obj);

            //Triangulo 2
            addPonto(x,-tamanho/2.0f,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, -1, 0, obj);
            addPonto(x,-tamanho/2.0f,z+count,obj);
            addTexPonto(u, v + v_step, obj);
            addNormPonto(0, -1, 0, obj);
            addPonto(x+count,-tamanho/2.0f,z+count,obj);
            addTexPonto(u - u_step, v + v_step, obj);
            addNormPonto(0, -1, 0, obj);
            v += v_step;
        }
        u -= u_step;
    }

    u = 0.0;
    v = 0.0;

    // Base frontal
    for (float x = -tamanho/2.0f; x < tamanho/2.0f; x+=count ) {
        v = 0.0;
        for(float z = -tamanho/2.0f; z < tamanho/2.0f; z+=count) {
            //Triangulo 1
            addPonto(x,tamanho/2.0f,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 1, 0, obj);
            addPonto(x+count,tamanho/2.0f,z,obj);
            addTexPonto(u + u_step, v, obj);
            addNormPonto(0, 1, 0, obj);
            addPonto(x+count,tamanho/2.0f,z+count,obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(0, 1, 0, obj);

            //Triangulo 2
            addPonto(x,tamanho/2.0f,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(0, 1, 0, obj);
            addPonto(x+count,tamanho/2.0f,z+count,obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(0, 1, 0, obj);
            addPonto(x,tamanho/2.0f,z+count,obj);
            addTexPonto(u, v + v_step, obj);
            addNormPonto(0, 1, 0, obj);
            v += v_step;
        }
        u += u_step;
    }

    u = 0.0;
    v = 0.0;

    // Base lateral esquerda
    for (float z = -tamanho/2.0f; z < tamanho/2.0f; z+=count ) {
        u = 0.0;
        for(float y = -tamanho/2.0f; y < tamanho/2.0f; y+=count) {
            //Triangulo 1
            addPonto(-tamanho/2.0f,y,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(-1, 0, 0, obj);
            addPonto(-tamanho/2.0f,y+count,z+count,obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(-1, 0, 0, obj);
            addPonto(-tamanho/2.0f,y,z+count,obj);
            addTexPonto(u, v + v_step, obj);
            addNormPonto(-1, 0, 0, obj);

            //Triangulo 2
            addPonto(-tamanho/2.0f,y,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(-1, 0, 0, obj);
            addPonto(-tamanho/2.0f,y+count,z,obj);
            addTexPonto(u + u_step, v, obj);
            addNormPonto(-1, 0, 0, obj);
            addPonto(-tamanho/2.0f,y+count,z+count,obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(-1, 0, 0, obj);
            u += u_step;
        }
        v += v_step;
    }

    u = 0.0;
    v = 0.0;

    // Base lateral direita
    for (float z = -tamanho/2.0f; z < tamanho/2.0f; z+=count ) {
        u = 1.0;
        for(float y = -tamanho/2.0f; y < tamanho/2.0f; y+=count) {
            //Triangulo 1
            addPonto(tamanho/2.0f,y,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(1, 0, 0, obj);
            addPonto(tamanho/2.0f,y,z+count,obj);
            addTexPonto(u, v + v_step, obj);
            addNormPonto(1, 0, 0, obj);
            addPonto(tamanho/2.0f,y+count,z+count,obj);
            addTexPonto(u - u_step, v + v_step, obj);
            addNormPonto(1, 0, 0, obj);

            //Triangulo 2
            addPonto(tamanho/2.0f,y,z,obj);
            addTexPonto(u, v, obj);
            addNormPonto(1, 0, 0, obj);
            addPonto(tamanho/2.0f,y+count,z+count,obj);
            addTexPonto(u - u_step, v + v_step, obj);
            addNormPonto(1, 0, 0, obj);
            addPonto(tamanho/2.0f,y+count,z,obj);
            addTexPonto(u - u_step, v, obj);
            addNormPonto(1, 0, 0, obj);
            u -= u_step;
        }
        v += v_step;
    }

    index_vertices(obj);
    makeFile(fich, obj);
}


void coneFunc(float radius, float height, int slices, int stacks, std::string fich) {
    
    OBJETO obj = new objeto();
    float angulo = (float)(2*M_PI) / ((float)slices);
    float anguloInicial = 0;
    float alturaBase = 0;//-height / 2.0;
    float rAtual = (float)radius;
    float rSeguinte = 0;
    float hAtual = (float)alturaBase;
    float hSeguinte = 0;
    float u = 0.0;
    float v = 0.0;
    float u_step = 1.0f / float(slices);
    float v_step = 1.0f / float(stacks);


    //base
    while(abs(anguloInicial - 2 * M_PI) > 0.01) {
        addPonto(radius * cos(anguloInicial + angulo), radius * sin(anguloInicial + angulo), alturaBase, obj);
        addTexPonto(u + u_step, v, obj);
        addNormPonto(0, 0, -1, obj);
        addPonto(0, 0, alturaBase, obj);
        addTexPonto(u + u_step/2, v + v_step, obj);
        addNormPonto(0, 0, -1, obj);
        addPonto(radius * cos(anguloInicial), radius * sin(anguloInicial), alturaBase, obj);
        addTexPonto(u, v, obj);
        addNormPonto(0, 0, -1, obj);
        anguloInicial += angulo;
        u += u_step;
    }
    float side_angle = atan(radius / height);
    u = 0.0;
    v = 0.0;
    //altura
    anguloInicial = 0;
    for (int i = 0; i < stacks; i++) {
        u = 1.0;
        rSeguinte = rAtual - (float)radius / (float)stacks;
        hSeguinte = hAtual + (float)height / (float)stacks;
        while (abs(anguloInicial - 2 * M_PI) > 0.01) {
            addPonto(rAtual * cos(anguloInicial), rAtual * sin(anguloInicial), hAtual, obj);
            addTexPonto(u, v, obj);
            addNormPonto(cos(anguloInicial), sin(anguloInicial), sin(side_angle), obj);
            addPonto(rSeguinte * cos(anguloInicial), rSeguinte * sin(anguloInicial), hSeguinte, obj);
            addTexPonto(u, v + v_step, obj);
            addNormPonto(cos(anguloInicial), sin(anguloInicial), sin(side_angle), obj);
            addPonto(rSeguinte * cos(anguloInicial + angulo), rSeguinte * sin(anguloInicial + angulo), hSeguinte, obj);
            addTexPonto(u - u_step, v + v_step, obj);
            addNormPonto(cos(anguloInicial + angulo), sin(anguloInicial + angulo), sin(side_angle), obj);

            addPonto(rAtual * cos(anguloInicial), rAtual * sin(anguloInicial), hAtual, obj);
            addTexPonto(u, v, obj);
            addNormPonto(cos(anguloInicial), sin(anguloInicial), sin(side_angle), obj);
            addPonto(rSeguinte * cos(anguloInicial + angulo), rSeguinte * sin(anguloInicial + angulo), hSeguinte, obj);
            addTexPonto(u - u_step, v + v_step, obj);
            addNormPonto(cos(anguloInicial + angulo), sin(anguloInicial + angulo), sin(side_angle), obj);
            addPonto(rAtual * cos(anguloInicial + angulo), rAtual * sin(anguloInicial + angulo), hAtual, obj);
            addTexPonto(u - u_step, v, obj);
            addNormPonto(cos(anguloInicial + angulo), sin(anguloInicial + angulo), sin(side_angle), obj);

            anguloInicial += angulo;
            u -= u_step;
        }
        v += v_step;
        hAtual = hSeguinte;
        rAtual = rSeguinte;
        anguloInicial = 0;
    }
    index_vertices(obj);
    makeFile(fich, obj);
}


void sphereFunc(float radius, int slices, int stacks, std::string fich){
    OBJETO obj = new objeto();
    float lat = 0.0;
    float lon = 0.0;
    float lat_step = (float) M_PI / (float) stacks;
    float lon_step = (float) 2 * M_PI / (float) slices;
    float u = 1.0;
    float v = 1.0;
    float u_step = 1.0 / float(slices);
    float v_step = 1.0 / float(stacks);
    
    while(abs(lat - M_PI) > 0.01){
        float next_lat = lat + lat_step;
        u = 1.0;

        while(abs(lon - 2 * M_PI) > 0.01){
            float next_lon = lon + lon_step;
            // First triangle
            addPonto(radius * cos(lon) * sin(lat), radius * sin(lon) * sin(lat), radius * cos(lat), obj);
            addTexPonto(u, v, obj);
            addNormPonto(cos(lon) * sin(lat), sin(lon) * sin(lat), cos(lat), obj);
            addPonto(radius * cos(next_lon) * sin(lat), radius * sin(next_lon) * sin(lat), radius * cos(lat), obj);
            addTexPonto(u - u_step, v, obj);
            addNormPonto(cos(next_lon) * sin(lat), sin(next_lon) * sin(lat), cos(lat), obj);
            addPonto(radius * cos(lon) * sin(next_lat), radius * sin(lon) * sin(next_lat), radius * cos(next_lat), obj);
            addTexPonto(u, v - v_step, obj);
            addNormPonto(cos(lon) * sin(next_lat), sin(lon) * sin(next_lat), cos(next_lat), obj);

            // Second triangle
            addPonto(radius * cos(next_lon) * sin(lat), radius * sin(next_lon) * sin(lat), radius * cos(lat), obj);
            addTexPonto(u - u_step, v, obj);
            addNormPonto(cos(next_lon) * sin(lat), sin(next_lon) * sin(lat), cos(lat), obj);
            addPonto(radius * cos(next_lon) * sin(next_lat), radius * sin(next_lon) * sin(next_lat), radius * cos(next_lat), obj);
            addTexPonto(u - u_step, v - v_step, obj);
            addNormPonto(cos(next_lon) * sin(next_lat), sin(next_lon) * sin(next_lat), cos(next_lat), obj);
            addPonto(radius * cos(lon) * sin(next_lat), radius * sin(lon) * sin(next_lat), radius * cos(next_lat), obj);
            addTexPonto(u, v - v_step, obj);
            addNormPonto(cos(lon) * sin(next_lat), sin(lon) * sin(next_lat), cos(next_lat), obj);
            
            u -= u_step;

            lon = next_lon;
        }
        v -= v_step;

        lon = 0.0;
        lat = next_lat;
    }
    index_vertices(obj);
    makeFile(fich, obj);
}

void inverseSphereFunc(float radius, int slices, int stacks, std::string fich){
    OBJETO obj = new objeto();
    float lat = 0.0;
    float lon = 0.0;
    float lat_step = (float) M_PI / (float) stacks;
    float lon_step = (float) 2 * M_PI / (float) slices;
    float u = 1.0;
    float v = 1.0;
    float u_step = 1.0 / float(slices);
    float v_step = 1.0 / float(stacks);
    
    while(abs(lat - M_PI) > 0.01){
        float next_lat = lat + lat_step;
        u = 1.0;

        while(abs(lon - 2 * M_PI) > 0.01){
            float next_lon = lon + lon_step;
            // First triangle
            addPonto(radius * cos(lon) * sin(lat), radius * sin(lon) * sin(lat), radius * cos(lat), obj);
            addTexPonto(u, v, obj);
            addNormPonto(cos(lon) * sin(lat), sin(lon) * sin(lat), cos(lat), obj);
            addPonto(radius * cos(lon) * sin(next_lat), radius * sin(lon) * sin(next_lat), radius * cos(next_lat), obj);
            addTexPonto(u, v - v_step, obj);
            addNormPonto(cos(lon) * sin(next_lat), sin(lon) * sin(next_lat), cos(next_lat), obj);
            addPonto(radius * cos(next_lon) * sin(lat), radius * sin(next_lon) * sin(lat), radius * cos(lat), obj);
            addTexPonto(u - u_step, v, obj);
            addNormPonto(cos(next_lon) * sin(lat), sin(next_lon) * sin(lat), cos(lat), obj);

            // Second triangle
            addPonto(radius * cos(next_lon) * sin(lat), radius * sin(next_lon) * sin(lat), radius * cos(lat), obj);
            addTexPonto(u - u_step, v, obj);
            addNormPonto(cos(next_lon) * sin(lat), sin(next_lon) * sin(lat), cos(lat), obj);
            addPonto(radius * cos(lon) * sin(next_lat), radius * sin(lon) * sin(next_lat), radius * cos(next_lat), obj);
            addTexPonto(u, v - v_step, obj);
            addNormPonto(cos(lon) * sin(next_lat), sin(lon) * sin(next_lat), cos(next_lat), obj);
            addPonto(radius * cos(next_lon) * sin(next_lat), radius * sin(next_lon) * sin(next_lat), radius * cos(next_lat), obj);
            addTexPonto(u - u_step, v - v_step, obj);
            addNormPonto(cos(next_lon) * sin(next_lat), sin(next_lon) * sin(next_lat), cos(next_lat), obj);
            
            u -= u_step;

            lon = next_lon;
        }
        v -= v_step;

        lon = 0.0;
        lat = next_lat;
    }
    index_vertices(obj);
    makeFile(fich, obj);
}


void torusFunc(float innerRadius, float outterRadius, int slices, int rings, std::string fich){
    OBJETO obj = new objeto();
    std::list<ponto> first_circle;
    float angle = 0;
    float slice_angle_step = (2 * M_PI) / (float) slices;
    float ring_angle_step = (2 * M_PI) / (float) rings;
    float radius = (outterRadius - innerRadius) / 2.0f;
    float u = 0.0;
    float v = 0.0;
    float u_step = 1.0f / float(slices);
    float v_step = 1.0f / float(rings);
    std::vector<ponto> circ_norms;

    // draw the first circle
    while(abs(angle - 2 * M_PI) > 0.01){
        ponto p;
        p.x = cos(angle) * radius + (outterRadius - radius);
        p.y = cos(angle) * radius + (outterRadius - radius);
        p.z = sin(angle) * radius;
        first_circle.push_back(p);
        ponto n;
        n.x = cos(angle);
        n.y = cos(angle);
        n.z = sin(angle);
        circ_norms.push_back(n);
        angle += slice_angle_step;
    }
    int i = 0;
    angle = 0;
    while(abs(angle - 2 * M_PI) > 0.01){
        u = 0.0;
        i = 0;
        for (std::list<ponto>::iterator p = first_circle.begin(); p != first_circle.end(); ++p){
            std::advance(p, 1);
            ponto n = circ_norms[i], n_next;
            if(p == first_circle.end()){
                p = first_circle.begin();
                n_next = circ_norms[0];
            }else n_next = circ_norms[i+1];
            std::list<ponto>::iterator p_next = p;
            if(p == first_circle.begin()){
                p = first_circle.end();
            }
            std::advance(p, -1);

            // Triangle 1
            addPonto(sin(angle) * p->x, cos(angle) * p->y, p->z, obj);
            addTexPonto(u, v, obj);
            addNormPonto(sin(angle) * n.x, cos(angle) * n.y, n.z, obj);

            addPonto(sin(angle + ring_angle_step) * p->x, cos(angle + ring_angle_step) * p->y, p->z, obj);
            addTexPonto(u, v + v_step, obj);
            addNormPonto(sin(angle + ring_angle_step) * n.x, cos(angle + ring_angle_step) * n.y, n.z, obj);
            
            addPonto(sin(angle + ring_angle_step) * p_next->x, cos(angle + ring_angle_step) * p_next->y, p_next->z, obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(sin(angle + ring_angle_step) * n_next.x, cos(angle + ring_angle_step) * n_next.y, n_next.z, obj);



            // Triangle 2
            addPonto(sin(angle) * p->x, cos(angle) * p->y, p->z, obj);
            addTexPonto(u, v, obj);
            addNormPonto(sin(angle) * n.x, cos(angle) * n.y, n.z, obj);
            
            addPonto(sin(angle + ring_angle_step) * p_next->x, cos(angle + ring_angle_step) * p_next->y, p_next->z, obj);
            addTexPonto(u + u_step, v + v_step, obj);
            addNormPonto(sin(angle + ring_angle_step) * n_next.x, cos(angle + ring_angle_step) * n_next.y, n_next.z, obj);
            
            addPonto(sin(angle) * p_next->x, cos(angle) * p_next->y, p_next->z, obj);
            addTexPonto(u + u_step, v, obj);
            addNormPonto(sin(angle) * n_next.x, cos(angle) * n_next.y, n_next.z, obj);
            
            u += u_step;
            i++;
        }
        v += v_step;
        angle += ring_angle_step;
    }
    index_vertices(obj);
    makeFile(fich, obj);
}


// Read patch file for bezier surfaces
void read_patch_file(std::string file_path, std::vector<PONTO> *vertices, std::vector<int> *indices){
    
    std::ifstream file(file_path);
	std::string text;

    getline(file, text);
	int surface_count = atoi(text.c_str());

    // Read indices
    for (int i = 0; i < surface_count; i++){
        getline(file, text);
        char* ind_str = strtok((char*)text.c_str(), " ,\n");
        indices->push_back(atoi(ind_str));
        for (int j = 0; j < 15; j++){
            ind_str = strtok(NULL, " ,\n");
            indices->push_back(atoi(ind_str));
        }
    }

    getline(file, text);
    int vertice_count = atoi(text.c_str());

    // Read vertices
    for (int i = 0; i < vertice_count; i++){
        getline(file, text);
        PONTO p = new ponto(); 
        char* vert_str = strtok((char*)text.c_str(), " ,\n");
        p->x = atof(vert_str);
        vert_str = strtok(NULL, " ,\n");
        p->y = atof(vert_str);
        vert_str = strtok(NULL, " ,\n");
        p->z = atof(vert_str);
        vertices->push_back(p);
    }
}


void getBezierPoint(float u, float v, float** matrixX, float** matrixY, float** matrixZ, float* pos) {
    float bezierMatrix[4][4] = { { -1.0f, 3.0f, -3.0f, 1.0f },
                               { 3.0f, -6.0f, 3.0f, 0.0f },
                               { -3.0f, 3.0f, 0.0f, 0.0f },
                               { 1.0f,  0.0f, 0.0f, 0.0f } };

    float vetorU[4] = { u * u * u, u * u, u, 1 };
    float vetorV[4] = { v * v * v, v * v, v, 1 };

    float vetorAux[4];
    float px[4];
    float py[4];
    float pz[4];

    float mx[4];
    float my[4];
    float mz[4];

    multMatrixVector((float*)bezierMatrix, vetorV, vetorAux);
    multMatrixVector((float*)matrixX, vetorAux, px);
    multMatrixVector((float*)matrixY, vetorAux, py);
    multMatrixVector((float*)matrixZ, vetorAux, pz);

    multMatrixVector((float*)bezierMatrix, px, mx);
    multMatrixVector((float*)bezierMatrix, py, my);
    multMatrixVector((float*)bezierMatrix, pz, mz);


    pos[0] = (vetorU[0] * mx[0]) + (vetorU[1] * mx[1]) + (vetorU[2] * mx[2]) + (vetorU[3] * mx[3]);
    pos[1] = (vetorU[0] * my[0]) + (vetorU[1] * my[1]) + (vetorU[2] * my[2]) + (vetorU[3] * my[3]);
    pos[2] = (vetorU[0] * mz[0]) + (vetorU[1] * mz[1]) + (vetorU[2] * mz[2]) + (vetorU[3] * mz[3]);

}


void getBezierNormalPoint(float u, float v, float** matrixX, float** matrixY, float** matrixZ, float* pos) {
    float bezierMatrix[4][4] = { { -1.0f, 3.0f, -3.0f, 1.0f },
                               { 3.0f, -6.0f, 3.0f, 0.0f },
                               { -3.0f, 3.0f, 0.0f, 0.0f },
                               { 1.0f,  0.0f, 0.0f, 0.0f } };

    float vetorU[4] = { u * u * u, u * u, u, 1 };
    float vetorV[4] = { v * v * v, v * v, v, 1 };

    float vetorDerivU[4] = { 3 * u * u, 2 * u, 1, 0 };
    float vetorDerivV[4] = { 3 * v * v, 2 * v, 1, 0 };

    float vetorAux[4];
    float px[4];
    float py[4];
    float pz[4];

    float mxAux[4];
    float myAux[4];
    float mzAux[4];

    float uRes[3];

    //Calcular u
    multMatrixVector((float*)bezierMatrix, vetorDerivV, vetorAux);
    multMatrixVector((float*)matrixX, vetorAux, px);
    multMatrixVector((float*)matrixY, vetorAux, py);
    multMatrixVector((float*)matrixZ, vetorAux, pz);

    multMatrixVector((float*)bezierMatrix, px, mxAux);
    multMatrixVector((float*)bezierMatrix, py, myAux);
    multMatrixVector((float*)bezierMatrix, pz, mzAux);

    uRes[0] = (mxAux[0] * vetorU[0]) + (mxAux[1] * vetorU[1]) + (mxAux[2] * vetorU[2]) + (mxAux[3] * vetorU[3]);
    uRes[2] = (mzAux[0] * vetorU[0]) + (mzAux[1] * vetorU[1]) + (mzAux[2] * vetorU[2]) + (mzAux[3] * vetorU[3]);
    uRes[1] = (myAux[0] * vetorU[0]) + (myAux[1] * vetorU[1]) + (myAux[2] * vetorU[2]) + (myAux[3] * vetorU[3]);


    //Calcular v
    float vetorAux2[4];
    float px2[4];
    float py2[4];
    float pz2[4];

    float mxAux2[4];
    float myAux2[4];
    float mzAux2[4];

    float vRes[3];

    multMatrixVector((float*)bezierMatrix, vetorV, vetorAux2);
    multMatrixVector((float*)matrixX, vetorAux2, px2);
    multMatrixVector((float*)matrixY, vetorAux2, py2);
    multMatrixVector((float*)matrixZ, vetorAux2, pz2);

    multMatrixVector((float*)bezierMatrix, px2, mxAux2);
    multMatrixVector((float*)bezierMatrix, py2, myAux2);
    multMatrixVector((float*)bezierMatrix, pz2, mzAux2);
    
    vRes[0] = (mxAux2[0] * vetorDerivU[0]) + (mxAux2[1] * vetorDerivU[1]) + (mxAux2[2] * vetorDerivU[2]) + (mxAux2[3] * vetorDerivU[3]);
    vRes[1] = (myAux2[0] * vetorDerivU[0]) + (myAux2[1] * vetorDerivU[1]) + (myAux2[2] * vetorDerivU[2]) + (myAux2[3] * vetorDerivU[3]);
    vRes[2] = (mzAux2[0] * vetorDerivU[0]) + (mzAux2[1] * vetorDerivU[1]) + (mzAux2[2] * vetorDerivU[2]) + (mzAux2[3] * vetorDerivU[3]);

    cross(uRes, vRes, pos);
    normalize(pos);

}



// Make bezier model
void make_bezier_model(std::vector<PONTO>* vertices, std::vector<int>* indices, int tessellation, std::string file){

    OBJETO obj = new objeto();
    
    float pos[4][3];
    float matrixX[4][4];
    float matrixY[4][4];
    float matrixZ[4][4];

    float u = 0.0;
    float v = 0.0;
    float step = 1 / (float)tessellation;

    for (size_t p = 0; p < (*indices).size(); p += 16) {
        for (size_t i = 0; i < tessellation; i++) {
            for (size_t j = 0; j < tessellation; j++) {

                for (size_t a = 0; a < 4; a++) {
                    for (size_t b = 0; b < 4; b++) {

                        matrixX[a][b] = (*vertices).at((*indices).at(p + a * 4 + b))->x;
                        matrixY[a][b] = (*vertices).at((*indices).at(p + a * 4 + b))->y;
                        matrixZ[a][b] = (*vertices).at((*indices).at(p + a * 4 + b))->z;

                    }
                }

                getBezierPoint(u, v, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[0]);
                getBezierPoint(u, v + step, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[1]);
                getBezierPoint(u + step, v, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[2]);
                getBezierPoint(u + step, v + step, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[3]);
                
                // Triangle 1
                addPonto(pos[3][0], pos[3][2], pos[3][1], obj);
                addTexPonto(u + step, v + step, obj);
                addPonto(pos[2][0], pos[2][2], pos[2][1], obj);
                addTexPonto(u + step, v, obj);
                addPonto(pos[0][0], pos[0][2], pos[0][1], obj);
                addTexPonto(u, v, obj);

                // Triangle 2
                addPonto(pos[0][0], pos[0][2], pos[0][1], obj);
                addTexPonto(u, v, obj);
                addPonto(pos[1][0], pos[1][2], pos[1][1], obj);
                addTexPonto(u, v + step, obj);
                addPonto(pos[3][0], pos[3][2], pos[3][1], obj);
                addTexPonto(u + step, v + step, obj);

                // Make normals
                getBezierNormalPoint(u, v, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[0]);
                getBezierNormalPoint(u, v + step, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[1]);
                getBezierNormalPoint(u + step, v, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[2]);
                getBezierNormalPoint(u + step, v + step, (float**)matrixX, (float**)matrixY, (float**)matrixZ, pos[3]);
    
                addNormPonto(pos[3][0], pos[3][2], pos[3][1], obj);
                addNormPonto(pos[2][0], pos[2][2], pos[2][1], obj);
                addNormPonto(pos[0][0], pos[0][2], pos[0][1], obj);
                addNormPonto(pos[0][0], pos[0][2], pos[0][1], obj);
                addNormPonto(pos[1][0], pos[1][2], pos[1][1], obj);
                addNormPonto(pos[3][0], pos[3][2], pos[3][1], obj);
                

                v += step;
            }
            v = 0.0;
            u += step;
        }
        v = 0.0;
        u = 0.0;
    }
    index_vertices(obj);
    makeFile(file, obj);
}


// Função main que recebe os argumentos passados no terminal
int main(int argc, char **argv) {

    if(argc == 0) {
        printf("Arguments invalid !!!");
        return 0;
    }

    if (strcmp(argv[1], "plane") == 0){
        if(argc != 5) {
            printf("Wrong number of arguments !!! ");
            return 0;
        }
        planeFunc(atof(argv[2]),atoi(argv[3]),argv[4]);
    }
    else if (strcmp(argv[1], "box") == 0){
        if(argc != 5) {
            printf("Wrong number of arguments !!! ");
            return 0;
        }
        boxFunc(atof(argv[2]),atoi(argv[3]),argv[4]);
    }
    else if (strcmp(argv[1], "cone") == 0){
        if(argc != 7) {
            printf("Wrong number of arguments !!! ");
            return 0;
        }
        coneFunc(atof(argv[2]),atof(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6]);
    }
    else if (strcmp(argv[1], "sphere") == 0){
        if(argc != 6) {
            printf("Wrong number of arguments !!! ");
            return 0;
        }
        sphereFunc(atof(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5]);
    }
    else if (strcmp(argv[1], "inverse_sphere") == 0){
        if(argc != 6) {
            printf("Wrong number of arguments !!! ");
            return 0;
        }
        inverseSphereFunc(atof(argv[2]),atoi(argv[3]),atoi(argv[4]),argv[5]);
    }
    else if (strcmp(argv[1], "torus") == 0){
        if(argc != 7) {
            printf("Wrong number of arguments !!! ");
            return 0;
        }
        torusFunc(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),argv[6]);
    }
    else if (strcmp(argv[1], "bezier") == 0 || strcmp(argv[1], "patch") == 0){
        if(argc != 5){
            printf("Wrong number of arguments !!!");
            return 0;
        }
        std::vector<PONTO> vertices;
        std::vector<int> indices;
        read_patch_file(argv[2], &vertices, &indices);
        make_bezier_model(&vertices, &indices, atoi(argv[3]), argv[4]);
    }

    return 0;
}