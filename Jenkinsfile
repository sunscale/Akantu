pipeline {
  agent {
    docker {
      image 'docker-slave'
    }
    
  }
  stages {
    stage('Configure') {
      agent {
        docker {
          image 'docker-slave'
        }
        
      }
      steps {
        sh 'makdir build && cd build'
        sh 'cmake ..'
      }
    }
  }
}