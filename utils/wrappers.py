import time

def logger( function ):

	def inner( *args, **kwargs ):

		before = time.time()

		result = function( *args, **kwargs )

		after = time.time()

		print "\tPerformed {func_name} over {length} values in {time} s\n".format(
			func_name = function.func_name,
			length = len( result ),
			time = after - before )

		return result

	return inner