
process say_hi {
  
  input:
  val x

  output:
  stdout

  script:
  """
  python ${baseDir}/modules/test_py/x1_test.py -c ${x}
  """
}