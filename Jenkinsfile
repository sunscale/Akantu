pipeline {
  parameters {string(defaultValue: '', description: 'api-token', name: 'API_TOKEN')
              string(defaultValue: '', description: 'buildable phid', name: 'TARGET_PHID')
              string(defaultValue: '', description: 'Commit id', name: 'COMMIT_ID')
              string(defaultValue: '', description: 'Diff id', name: 'DIFF_ID')
  }

  options {
    disableConcurrentBuilds()
  }

  agent {
    dockerfile {
      additionalBuildArgs '--tag akantu-environment'
    }
  }
  stages {
    stage('Configure') {
      steps {
        sh """
        env
        mkdir -p build
        cd build
        cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE \
              -DAKANTU_IMPLICIT:BOOL=TRUE \
              -DAKANTU_PARALLEL:BOOL=TRUE \
              -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE \
              -DAKANTU_TESTS:BOOL=TRUE ..
        """
      }
	    post {
				failure {
					deleteDir()
				}
			}
    }
    stage('Compile') {
      steps {
				sh 'make -C build/src || true'
      }
    }

    stage ('Warnings gcc') {
      steps {
        warnings(consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']])
      }
    }

    stage('Compile python') {
      steps {
        sh 'make -C build/python || true'
      }
    }

    stage('Compile tests') {
      steps {
        sh 'make -C build/test || true'
      }
    }

    stage('Tests') {
      steps {
        sh '''
          rm -rf build/gtest_reports
          cd build/
          #source ./akantu_environement.sh
        
          ctest -T test --no-compress-output || true
        '''
				
        def TAG = sh returnStdout: true, script: 'head -n 1 < build/Testing/TAG'

				if (fileExists("build/Testing/${TAG}/Test.xml")) {
					sh "cp build/Testing/${TAG}/Test.xml CTestResults.xml"
        }
      }
    }
  }
  environment {
    BLA_VENDOR = 'OpenBLAS'
    OMPI_MCA_plm = 'isolated'
    OMPI_MCA_btl = 'tcp,self'
  }
  post {
    always {
			step([$class: 'XUnitBuilder',
			  thresholds: [
          [$class: 'SkippedThreshold', failureThreshold: '0'],
          [$class: 'FailedThreshold', failureThreshold: '0']],
				tools: [
			  	[$class: 'CTestType', pattern: 'CTestResults.xml', skipNoTestFiles: true]
				]])
      step([$class: 'XUnitBuilder',
         thresholds: [
          [$class: 'SkippedThreshold', failureThreshold: '100'],
          [$class: 'FailedThreshold', failureThreshold: '0']],
         tools: [
					[$class: 'GoogleTestType', pattern: 'build/gtest_reports/**', skipNoTestFiles: true]
				]])
      archiveArtifacts artifacts: 'build/Testing/**', fingerprint: true
      createArtifact()
    }

    success {
      sendFailPass('pass')
    }

    failure {
      emailext(
        body: '''${SCRIPT, template="groovy-html.template"}''',
				mimeType: 'text/html',
        subject: "[Jenkins] ${currentBuild.fullDisplayName} Failed",
				recipientProviders: [[$class: 'CulpritsRecipientProvider']],
				to: 'akantu-admins@akantu.ch',
				replyTo: 'akantu-admins@akantu.ch',
				attachLog: true,
        compressLog: false)
      sendFailPass('fail')
    }
  }
}

def sendFailPass(state) {
    sh """
       set +x
       curl https://c4science.ch/api/harbormaster.sendmessage \
            -d api.token=${API_TOKEN} \
            -d buildTargetPHID=${TARGET_PHID} \
            -d type=${state}
       """
}

def createArtifact() {
    sh """ set +x
       curl https://c4science.ch/api/harbormaster.createartifact \
            -d api.token=${API_TOKEN} \
            -d buildTargetPHID=${TARGET_PHID} \
            -d artifactKey="Jenkins URI" \
            -d artifactType=uri \
            -d artifactData[uri]=${BUILD_URL} \
            -d artifactData[name]="View Jenkins result" \
            -d artifactData[ui.external]=1
       """
}
