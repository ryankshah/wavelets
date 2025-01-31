# lib/wavelets/ndwt.ex
defmodule Wavelets.NDWT do
  @moduledoc """
  Implementation of n-dimensional Discrete Wavelet Transform
  """

  alias Wavelets.{DWT, Filter}

  @doc """
  Performs n-dimensional forward discrete wavelet transform.
  The input signal should be a nested list structure matching the dimensions.
  Returns a tuple with approximation and a list of detail coefficients.
  """
  @spec forward(list(), Filter.t(), non_neg_integer(), Wavelets.precision()) ::
          {list(), list()}
  def forward(signal, filter, dims, _precision \\ :double) do
    case dims do
      1 ->
        DWT.forward_1d(signal, filter)

      2 ->
        DWT.forward_2d(signal, filter)

      _ ->
        do_forward_nd(signal, filter, dims)
    end
  end

  @doc """
  Performs n-dimensional inverse discrete wavelet transform
  """
  @spec inverse(
          list(),
          list(),
          Filter.t(),
          non_neg_integer(),
          Wavelets.precision()
        ) :: list()
  def inverse(approximation, details, filter, dims, _precision \\ :double) do
    case dims do
      1 ->
        Wavelets.IDWT.inverse_1d(approximation, details, filter)

      2 ->
        Wavelets.IDWT.inverse_2d(approximation, details, filter)

      _ ->
        do_inverse_nd(approximation, details, filter, dims)
    end
  end

  # Private helper functions
  defp do_forward_nd(signal, filter, dims) do
    # Helper function to apply transform along a specific dimension
    transform_dim = fn signal_nd, dim ->
      case dim do
        0 ->
          # Base case: apply 1D transform to the innermost lists
          Enum.map(signal_nd, &DWT.forward_1d(&1, filter))

        _ ->
          # Recursive case: traverse to inner dimensions
          Enum.map(signal_nd, &transform_dim.(&1, dim - 1))
      end
    end

    # Apply transforms along each dimension
    coeffs =
      Enum.reduce((dims - 1)..0, signal, fn dim, acc ->
        transform_dim.(acc, dim)
      end)

    # Separate approximation and details
    extract_coeffs(coeffs, dims)
  end

  defp do_inverse_nd(approximation, details, filter, dims) do
    # Reconstruct the coefficient structure
    coeffs = merge_coeffs(approximation, details, dims)

    # Helper function to apply inverse transform along a specific dimension
    inverse_dim = fn signal_nd, dim ->
      case dim do
        0 ->
          # Base case: apply 1D inverse transform
          Enum.map(signal_nd, fn {a, d} ->
            Wavelets.IDWT.inverse_1d(a, d, filter)
          end)

        _ ->
          # Recursive case: traverse to inner dimensions
          Enum.map(signal_nd, &inverse_dim.(&1, dim - 1))
      end
    end

    # Apply inverse transforms along each dimension
    Enum.reduce(0..(dims - 1), coeffs, fn dim, acc ->
      inverse_dim.(acc, dim)
    end)
  end

  defp extract_coeffs(coeffs, dims) do
    case dims do
      1 ->
        Enum.unzip(coeffs)

      _ ->
        {approx, details} =
          Enum.map(coeffs, &extract_coeffs(&1, dims - 1))
          |> Enum.unzip()

        {approx, List.flatten(details)}
    end
  end

  defp merge_coeffs(approximation, details, dims) do
    case dims do
      1 ->
        Enum.zip(approximation, details)

      _ ->
        detail_chunks = Enum.chunk_every(details, div(length(details), dims))

        Enum.zip(approximation, detail_chunks)
        |> Enum.map(fn {a, d} -> merge_coeffs(a, d, dims - 1) end)
    end
  end
end
