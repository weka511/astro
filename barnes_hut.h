Node* add(Body* body, Node* node);

void accelerationOn( Body const* body, Node const* node, double theta,
                     double& ax, double& ay, double& az);

void verlet( std::vector<Body*>& bodies, Node* root,
             double theta, double G, double dt );

void save_bodies( std::vector<Body*>& bodies, int i, std::string path);
			 