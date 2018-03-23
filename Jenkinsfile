pipeline {
  agent {
    docker {
      image 'richart/akantu-public-debian'
    }
    
  }
  stages {
    stage('Configure') {
      agent {
        docker {
          image 'richart/akantu-public-debian'
        }
        
      }
      steps {
        sh 'makdir build && cd build'
        sh 'cmake ..'
      }
    }
  }
}