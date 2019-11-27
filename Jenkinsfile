pipeline {
  parameters {string(defaultValue: '', description: 'api-token', name: 'API_TOKEN')
              string(defaultValue: '', description: 'buildable phid', name: 'BUILD_TARGET_PHID')
              string(defaultValue: '', description: 'Commit id', name: 'COMMIT_ID')
              string(defaultValue: '', description: 'Diff id', name: 'DIFF_ID')
	      string(defaultValue: 'PHID-PROJ-5eqyu6ooyjktagbhf473', description: 'ID of the project', name: 'PROJECT_ID')
  }

  options {
    disableConcurrentBuilds()
    //skipDefaultCheckout(true)
  }

  environment {
    PHABRICATOR_HOST = 'https://c4science.ch/api/'
    PYTHONPATH = sh returnStdout: true, script: 'echo ${WORKSPACE}/test/ci/script/'
    BLA_VENDOR = 'OpenBLAS'
    OMPI_MCA_plm = 'isolated'
    OMPI_MCA_btl = 'tcp,self'
  }
  
  agent {
    dockerfile {
      filename 'Dockerfile'
      dir 'test/ci/debian.testing'
      additionalBuildArgs '--tag akantu-environment'
    }
  }
  
  stages {
    stage('Checkout proper commit') {
      steps {
        checkout scm:  [$class: 'GitSCM',
          branches: [[name: "${COMMIT_ID}" ]]
        ], changelog: true
      }
    }
        
    stage('Lint') {
      steps {
        sh """
           arc lint --output json --rev HEAD^ | jq . -srM | tee lint.json
           ./test/ci/scripts/hbm send-arc-lint -f lint.json
           """
      }
    }
    
    stage('Configure') {
      steps {
        sh """#!/bin/bash
           set -o pipefail
           mkdir -p build
           cd build
           cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE \
                 -DAKANTU_IMPLICIT:BOOL=TRUE \
                 -DAKANTU_PARALLEL:BOOL=TRUE \
                 -DAKANTU_STRUCTURAL_MECHANICS:BOOL=TRUE \
                 -DAKANTU_HEAT_TRANSFER:BOOL=TRUE \
                 -DAKANTU_DAMAGE_NON_LOCAL:BOOL=TRUE \
                 -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE \
                 -DAKANTU_EXAMPLES:BOOL=TRUE \
                 -DAKANTU_BUILD_ALL_EXAMPLES:BOOL=TRUE \
                 -DAKANTU_TEST_EXAMPLES:BOOL=FALSE \
                 -DAKANTU_TESTS:BOOL=TRUE .. 2>&1 | tee configure.txt
           """
      }
      post {
        failure {
          uploadArtifact('configure.txt', 'Configure')
          deleteDir()
        }
      }
    }
    
    stage('Compile') {
      steps {
        sh '''#!/bin/bash
           set -o pipefail
           make -C build/src | tee compilation.txt
           '''
      }
      post {
        failure {
          uploadArtifact('compilation.txt', 'Compilation')
        }
      }
    }

    stage ('Warnings gcc') {
      steps {
        warnings(consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']])
      }
    }

    stage('Compile python') {
      steps {
        sh '''#!/bin/bash
           set -o pipefail

           make -C build/python | tee compilation_python.txt
           '''
      }
      post {
        failure {
          uploadArtifact('compilation_python.txt', 'Compilation_Python')
        }
      }
    }

    stage('Compile tests') {
      steps {
        sh '''#!/bin/bash
           set -o pipefail

           make -C build/test | tee compilation_test.txt
           '''
      }
      post {
        failure {
          uploadArtifact('compilation_test.txt', 'Compilation_Tests')
        }
      }
    }

    stage('Tests') {
      steps {
        sh '''
          #rm -rf build/gtest_reports
          cd build/
          #source ./akantu_environement.sh
        
          ctest -T test --no-compress-output || true
          tag=$(head -n 1 < Testing/TAG)
          if [ -e Testing/${tag}/Test.xml ]; then
            cp Testing/${tag}/Test.xml ../CTestResults.xml
          fi
        '''
      }
      //post {
      //  failure {
      //    zip zipFile: 'build.zip',  dir: 'build/', archive: true
      //  }
      //}
    }
  }
  post {
    always {
      createArtifact("./CTestResults.xml")

      step([$class: 'XUnitBuilder',
      thresholds: [
          [$class: 'SkippedThreshold', failureThreshold: '0'],
          [$class: 'FailedThreshold', failureThreshold: '0']],
      tools: [
        [$class: 'CTestType', pattern: 'CTestResults.xml', skipNoTestFiles: true]
      ]])
    }

    success {
      passed()
    }

    failure {
      failed()
    }
  }
}

def failed() {
  sh "./test/ci/scripts/hbm failed"
}

def passed() {
  sh "./test/ci/scripts/hbm passed"
}

def createArtifact(filename) {
  sh "./test/ci/scripts/hbm send-uri -k 'Jenkins URI' -u ${BUILD_URL} -l 'View Jenkins result'"
  sh "./test/ci/scripts/hbm send-ctest-results -f ${filename}"
}

def uploadArtifact(artifact, name) {
  sh "./test/ci/scripts/hbm upload-file -f ${artifact} -n \"${name}\" -v ${PROJECT_ID}"
}
